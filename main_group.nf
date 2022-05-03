#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["atlas_config":"$params.atlas_config",
                "atlas_directory":"$params.atlas_directory",
                "atlas_centroids":"$params.atlas_centroids",
                "run_average_bundles":"$params.run_average_bundles",
                "multi_parameters":"$params.multi_parameters",
                "minimal_vote_ratio":"$params.minimal_vote_ratio",
                "wb_clustering_thr":"$params.wb_clustering_thr",
                "seeds":"$params.seeds",
                "outlier_alpha":"$params.outlier_alpha",
                "register_processes":"$params.register_processes",
                "rbx_processes":"$params.rbx_processes",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "SCIL RecobundlesX pipeline"
log.info "==========================="
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
log.info "Atlas Config: $params.atlas_config"
log.info "Atlas Anat: $params.atlas_anat"
log.info "Atlas Directory: $params.atlas_directory"
log.info "Atlas Centroids: $params.atlas_centroids"
log.info ""
log.info "[Recobundles options]"
log.info "Multi-Parameters Executions: $params.multi_parameters"
log.info "Minimal Vote Percentage: $params.minimal_vote_ratio"
log.info "Whole Brain Clustering Threshold: $params.wb_clustering_thr"
log.info "Random Seeds: $params.seeds"
log.info "Outlier Removal Alpha: $params.outlier_alpha"
log.info ""
log.info ""

log.info "Input: $params.input"
root = file(params.input)
/* Watch out, files are ordered alphabetically in channel */
tractogram_for_recognition_group = Channel
     .fromFilePairs("$root/**/{*tracking*.*,}",
                    size: -1,
                    maxDepth:1) {it.parent.name}

Channel
    .fromPath("$root/**/*fa.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{anat_for_registration;anat_for_reference_centroids;anat_for_reference_bundles_group;anat_for_reference_bundles}

Channel
    .fromPath("$params.atlas_directory/*.nii.gz",
                    maxDepth:1)
    .into{atlas_anat_for_registration;atlas_anat_for_average}

Channel
    .fromPath("$params.atlas_directory/*_group.json",
                    maxDepth:1)
    .map{it -> [it.baseName, it]}
    .set{atlas_configs}

atlas_config_group = Channel.fromPath("$params.atlas_directory/config_groups.json")
Channel
    .fromPath("$params.atlas_directory/atlas/",
                    maxDepth:1)
    .into{atlas_directory_group;atlas_directory}

if (params.atlas_centroids) {
    atlas_centroids = Channel.fromPath("$params.atlas_centroids/*_centroid.trk")
}
else {
    atlas_centroids = Channel.empty()
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

anat_for_registration
    .combine(atlas_anat_for_registration)
    .set{anats_for_registration}
process Register_Anat {
    cpus params.register_processes
    memory '2 GB'

    input:
    set sid, file(native_anat), file(atlas) from anats_for_registration

    output:
    set sid, "${sid}__output0GenericAffine.mat" into transformation_for_recognition_group, transformation_for_recognition, transformation_for_centroids
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
    .set{anat_and_transformation}
process Transform_Centroids {
    input:
    set sid, file(anat), file(transfo) from anat_and_transformation
    each file(centroid) from atlas_centroids
    output:
    file "${sid}__${centroid.baseName}.trk"
    script:
    """
    scil_apply_transform_to_tractogram.py ${centroid} ${anat} ${transfo} ${sid}__${centroid.baseName}.trk --inverse --cut_invalid
    """
}


tractogram_for_recognition_group
    .join(anat_for_reference_bundles_group)
    .join(transformation_for_recognition_group)
    .combine(atlas_config_group)
    .combine(atlas_directory_group)
    .set{tractogram_and_transformation_group}
process Recognize_Bundles_Group {
    cpus params.rbx_processes
    memory params.rbx_memory_limit

    input:
    set sid, file(tractograms), file(reference), file(transfo), file(config), file(directory) from tractogram_and_transformation_group
    output:
    set sid, "*.trk" into bundles_for_sub_reco, tmp
    file "results.json"
    file "logfile.txt"
    script:
    """
    if [ `echo $tractograms | wc -w` -gt 1 ]; then
        scil_streamlines_math.py lazy_concatenate $tractograms tracking_concat.trk
    else
        mv $tractograms tracking_concat.trk
    fi
    scil_remove_invalid_streamlines.py tracking_concat.trk tractogram_ic.trk --reference ${reference} --remove_single_point --remove_overlapping_points
    mkdir tmp/
    scil_recognize_multi_bundles.py tractogram_ic.trk ${config} ${directory}/*/ ${transfo} --inverse --out_dir tmp/ \
        --log_level DEBUG --multi_parameters 8 --minimal_vote_ratio 0.5 \
        --tractogram_clustering_thr 12 15 18 --seeds 0 --processes $params.rbx_processes
    rm tractogram_ic.trk tracking_concat.trk
    mv tmp/* ./
    """
}


bundles_for_sub_reco
    .transpose()
    .combine(anat_for_reference_bundles, by: 0)
    .combine(transformation_for_recognition, by: 0)
    .combine(atlas_directory)
    .map{it -> [it[1].baseName, it[0], it[1], it[2], it[3], it[4]]}
    .combine(atlas_configs, by:0)
    .set{tractogram_and_transformation_config}
process Recognize_Bundles {
    cpus params.rbx_processes
    memory params.rbx_memory_limit

    input:
    set gName, sid, file(tractogram), file(reference), file(transfo), file(directory), file(config) from tractogram_and_transformation_config
    output:
    set sid, "*.trk" optional true into bundles_for_cleaning
    file "results.json" optional true 
    file "logfile.txt" optional true 
    script:
    """
    mkdir tmp/
    scil_recognize_multi_bundles.py ${tractogram} ${config} ${directory}/*/ ${transfo} --inverse --out_dir tmp/ \
        --log_level DEBUG --multi_parameters $params.multi_parameters --minimal_vote_ratio $params.minimal_vote_ratio \
        --tractogram_clustering_thr $params.wb_clustering_thr --seeds $params.seeds --processes $params.rbx_processes
    mv tmp/* ./
    """
}

bundles_for_cleaning
    .transpose()
    .set{all_bundles_for_cleaning}
process Clean_Bundles {
    input:
    set sid, file(bundle) from all_bundles_for_cleaning
    output:
    set sid, val(bname), "${sid}__*_cleaned.trk" optional true into bundle_for_density
    script:
    bname = bundle.name.take(bundle.name.lastIndexOf('.'))
    """
    scil_outlier_rejection.py ${bundle} "${sid}__${bname}_cleaned.trk" --alpha $params.outlier_alpha
    """
}

bundle_for_density
    .combine(transformation_for_average, by:0)
    .combine(atlas_anat_for_average)
    .set{all_bundles_transfo_for_average}
process Compute_Density_Bundles {
    input:
    set sid, val(bname), file(bundle), file(transfo), file(atlas) from all_bundles_transfo_for_average
    output:
    set bname, "*.nii.gz" into bundle_for_average
    when:
    params.run_average_bundles
    script:
    """
    scil_apply_transform_to_tractogram.py $bundle $atlas $transfo tmp.trk --remove_invalid
    scil_compute_streamlines_density_map.py tmp.trk "${sid}__${bname}_density.nii.gz"
    scil_image_math.py lower_threshold "${sid}__${bname}_density.nii.gz" 1 "${sid}__${bname}_binary.nii.gz"
    """
}


bundle_for_average
    .flatMap{ sid, bundles -> bundles.collect{[sid, it]} }
    .groupTuple(by: 0)
    .set{all_bundle_for_average}
process Average_Bundles {
    publishDir = params.Average_Bundles_Publish_Dir
    input:
    set val(bname), file(bundles_bin) from all_bundle_for_average
    output:
    file "${bname}_density.nii.gz"
    file "${bname}_binary.nii.gz"
    script:
    """
    scil_image_math.py addition *_density.nii.gz 0 ${bname}_density.nii.gz
    scil_image_math.py addition *_binary.nii.gz 0 ${bname}_binary.nii.gz
    """
}
