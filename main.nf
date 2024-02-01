#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["atlas_directory":"$params.atlas_directory",
                "run_average_bundles":"$params.run_average_bundles",
                "minimal_vote_ratio":"$params.minimal_vote_ratio",
                "seed":"$params.seed",
                "outlier_alpha":"$params.outlier_alpha",
                "register_processes":"$params.register_processes",
                "rbx_processes":"$params.rbx_processes",
                "single_dataset_size_GB":"$params.single_dataset_size_GB",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "SCIL RecobundlesX pipeline"
log.info "=========================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[Atlas]"
log.info "Atlas Directory: $params.atlas_directory"
log.info ""
log.info "[Recobundles options]"
log.info "Minimal Vote Percentage: $params.minimal_vote_ratio"
log.info "Random Seed: $params.seed"
log.info "Outlier Removal Alpha: $params.outlier_alpha"
log.info "Run Average Bundles: $params.run_average_bundles"
log.info ""
log.info ""

log.info "Input: $params.input"
root = file(params.input)
/* Watch out, files are ordered alphabetically in channel */
tractogram_for_recognition = Channel
     .fromFilePairs("$root/**/{*tracking*.*,}",
                    size: -1,
                    maxDepth:1) {it.parent.name}

Channel
    .fromPath("$root/**/*fa.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{anat_for_registration;anat_for_reference_centroids;anat_for_reference_bundles}

atlas_directory = Channel.fromPath("$params.atlas_directory/atlas")

Channel.fromPath("$params.atlas_directory/mni_masked.nii.gz")
    .into{atlas_anat;atlas_anat_for_average}
atlas_config = Channel.fromPath("$params.atlas_directory/config_fss_1.json")
centroids_dir = Channel.fromPath("$params.atlas_directory/centroids/",
                                   type: 'dir')
centroids_dir
    .into{atlas_centroids;atlas_centroids_for_average}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

anat_for_registration
    .combine(atlas_anat)
    .set{anats_for_registration}
process Register_Anat {
    cpus params.register_processes
    memory '2 GB'

    input:
    set sid, file(native_anat), file(atlas) from anats_for_registration

    output:
    set sid, "${sid}__output0GenericAffine.mat" into transformation_for_recognition, transformation_for_centroids
    set sid, "${sid}__output0GenericAffine.mat" into transformation_for_average
    file "${sid}__outputWarped.nii.gz"
    file "${sid}__native_anat.nii.gz"

    script:
    """
    export ANTS_RANDOM_SEED=1234
    antsRegistrationSyNQuick.sh -d 3 -f ${native_anat} -m ${atlas} -n ${params.register_processes} -o ${sid}__output -t a
    cp ${native_anat} ${sid}__native_anat.nii.gz
    """
}

anat_for_reference_centroids
    .join(transformation_for_centroids, by: 0)
    .combine(atlas_centroids)
    .set{anat_and_transformation}
process Transform_Centroids {
    input:
    set sid, file(anat), file(transfo), file(centroids_dir) from anat_and_transformation

    output:
    file "${sid}__*.trk" optional true

    script:
    """
    for centroid in ${centroids_dir}/*.trk;
        do bname=\${centroid/_centroid/}
        bname=\$(basename \$bname .trk)

        scil_apply_transform_to_tractogram.py \${centroid} ${anat} ${transfo} tmp.trk --inverse --keep_invalid -f
        scil_remove_invalid_streamlines.py tmp.trk ${sid}__\${bname}.trk --cut_invalid --remove_single_point --remove_overlapping_points --no_empty
    done
    """
}


tractogram_for_recognition
    .join(anat_for_reference_bundles)
    .join(transformation_for_recognition)
    .combine(atlas_config)
    .combine(atlas_directory)
    .set{tractogram_and_transformation}
process Recognize_Bundles {
    cpus params.rbx_processes
    memory { params.single_dataset_size_GB.GB * params.rbx_processes }

    input:
    set sid, file(tractograms), file(refenrence), file(transfo), file(config), file(directory) from tractogram_and_transformation

    output:
    set sid, "*.trk" into bundles_for_cleaning
    file "results.json"
    file "logfile.txt"

    script:
    """
    mkdir tmp/
    scil_recognize_multi_bundles.py ${tractograms} ${config} ${directory}/ ${transfo} --inverse --out_dir tmp/ \
        --log_level DEBUG --minimal_vote_ratio $params.minimal_vote_ratio \
        --seed $params.seed --processes $params.rbx_processes
    mv tmp/* ./
    """
}


bundles_for_cleaning
    .combine(transformation_for_average, by:0)
    .combine(atlas_anat_for_average)
    .set{all_bundles_transfo_for_clean_average}
process Clean_Bundles {
    input:
    set sid, file(bundles), file(transfo), file(atlas) from all_bundles_transfo_for_clean_average

    output:
    set sid, "${sid}__*_cleaned.trk" 
    file "${sid}__*.nii.gz" optional true into bundle_for_average

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi

        scil_outlier_rejection.py \${bundle} "${sid}__\${bname}_cleaned.trk" \
            --alpha $params.outlier_alpha
            
        if [ -s "${sid}__\${bname}_cleaned.trk" ]; then 
            if ${params.run_average_bundles}; then
                scil_apply_transform_to_tractogram.py "${sid}__\${bname}_cleaned.trk" \
		            ${atlas} ${transfo} tmp.trk --remove_invalid -f
            
                scil_compute_streamlines_density_map.py tmp.trk "${sid}__\${bname}_density_mni.nii.gz"

                scil_image_math.py lower_threshold "${sid}__\${bname}_density_mni.nii.gz" 0.01 \
                    "${sid}__\${bname}_binary_mni.nii.gz"
            fi
        else
            echo "After cleaning \${bundle} all streamlines were outliers."
        fi
    done
    """
}


bundle_for_average
    .flatten()
    .transpose()
    .combine(atlas_centroids_for_average)
    .groupTuple(by: 1)
    .set{all_bundle_for_average}
process Average_Bundles {
    publishDir = params.Average_Bundles_Publish_Dir
    input:
    set file(bundles), file(centroids_dir) from all_bundle_for_average

    output:
    file "*_average_density_mni.nii.gz" optional true
    file "*_average_binary_mni.nii.gz" optional true

    when:
        params.run_average_bundles

    script:
    """
    shopt -s nullglob
    mkdir tmp/
    for centroid in $centroids_dir/*.trk;
        do bname=\${centroid/_centroid/}
        bname=\$(basename \$bname .trk)

        nfiles=\$(find ./ -maxdepth 1 -type f -name "*__\${bname}_density_mni.nii.gz" | wc -l)
        if [[ \$nfiles -gt 0 ]];
            then scil_image_math.py addition *__\${bname}_density_mni.nii.gz 0 tmp/\${bname}_average_density_mni.nii.gz
            scil_image_math.py addition *__\${bname}_binary_mni.nii.gz 0 tmp/\${bname}_average_binary_mni.nii.gz
        fi
    done
    mv tmp/* ./
    """
}
