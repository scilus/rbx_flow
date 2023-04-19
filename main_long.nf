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
log.info ""
log.info "[Dispatch options]"
log.info "Register processes: $params.register_processes"
log.info "RBx processes: $params.rbx_processes"
log.info "Single dataset size: $params.single_dataset_size_GB"
log.info ""
log.info ""

log.info "Input: $params.input"
root = file(params.input)
/* Watch out, files are ordered alphabetically in channel */
Channel
     .fromPath("$root/**/*/*tracking*.*",
                    maxDepth:2)
    .map{[it.parent.parent.name, it.parent.name, it]}
    .groupTuple(by: [0,1])
    .into{tractograms_for_recognition_first; tractograms_for_recognition_second}

Channel
     .fromPath("$root/**/*/*fa.nii.gz",
                    maxDepth:2)
    .map{[it.parent.parent.name, it.parent.name, it]}
    .into{anat_for_registration;anat_for_reference_centroids;
            anat_for_reference_bundles_first;anat_for_reference_bundles_second}

if (!(params.atlas_directory)) {
    error "You must specify --atlas_directory, it is mandatory."
}

atlas_directory = Channel.fromPath("$params.atlas_directory/atlas")
Channel.fromPath("$params.atlas_directory/mni_masked.nii.gz")
    .into{atlas_anat;atlas_anat_for_revert;atlas_anat_for_average}
atlas_config_first = Channel.fromPath("$params.atlas_directory/config_fss_1.json")
atlas_config_second = Channel.fromPath("$params.atlas_directory/config_fss_2.json")
atlas_centroids = Channel.fromPath("$params.atlas_directory/centroids/*_centroid.trk")

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

anat_for_registration
    .combine(atlas_anat)
    .set{anats_for_registration}
process Register_Anat {
    publishDir = {"./results_rbx/$sid/$sess/$task.process/"}
    cpus params.register_processes
    memory '2 GB'

    input:
    set sid, sess, file(native_anat), file(atlas) from anats_for_registration

    output:
    set sid, sess, "${sid}-${sess}__output0GenericAffine.mat" into transformation_for_recognition_first,
        transformation_for_recognition_second, transformation_for_centroids
    set sid, sess, "${sid}-${sess}__output0GenericAffine.mat" into transformation_for_average
    file "${sid}-${sess}__outputWarped.nii.gz"
    file "${sid}-${sess}__native_anat.nii.gz"
    script:
    """
    export ANTS_RANDOM_SEED=1234
    antsRegistrationSyNQuick.sh -d 3 -f ${native_anat} -m ${atlas} \
        -n ${params.register_processes} -o ${sid}-${sess}__output -t a
    cp ${native_anat} ${sid}-${sess}__native_anat.nii.gz
    """
}


anat_for_reference_centroids
    .join(transformation_for_centroids, by: [0,1])
    .set{anat_and_transformation}
process Transform_Centroids {
    publishDir = {"./results_rbx/$sid/$sess/$task.process/"}

    input:
    set sid, sess, file(anat), file(transfo) from anat_and_transformation
    each file(centroid) from atlas_centroids

    output:
    file "${sid}-${sess}__${centroid.baseName}.trk" optional true
    script:
    """
    scil_apply_transform_to_tractogram.py ${centroid} ${anat} ${transfo} tmp.trk \
        --inverse --keep_invalid
    scil_remove_invalid_streamlines.py tmp.trk ${sid}-${sess}__${centroid.baseName}.trk \
        --cut_invalid --remove_single_point --remove_overlapping_points --no_empty
    """
}


tractograms_for_recognition_first
    .join(anat_for_reference_bundles_first, by: [0,1])
    .join(transformation_for_recognition_first, by: [0,1])
    .combine(atlas_config_first)
    .combine(atlas_directory)
    .combine(atlas_anat_for_revert)
    .set{tractograms_and_transformation_atlas}
process Recognize_Bundles_Intial_MNI {
    publishDir = {"./results_rbx/$sid/$sess/$task.process/"}
    memory { params.single_dataset_size_GB.GB * params.rbx_processes }
    cpus params.rbx_processes

    input:
    set sid, sess, file(tractograms), file(reference), file(transfo),
        file(config), file(directory), file(atlas_anat) from tractograms_and_transformation_atlas
    
    output:
    set sid, sess, "${sid}-${sess}__recon_mni/" into bundles_for_longitudinal
    script:
    """
    mkdir tmp/ ${sid}-${sess}__recon_mni/
    scil_recognize_multi_bundles.py $tractograms ${config} ${directory}/ \
        ${transfo} --inverse --out_dir tmp/ \
        --log_level DEBUG --minimal_vote_ratio $params.minimal_vote_ratio \
        --seed $params.seed --processes $params.rbx_processes
    for i in tmp/*.trk;
        do bname=\$(basename \${i})
        scil_apply_transform_to_tractogram.py \${i} ${atlas_anat} \
        ${transfo} ${sid}-${sess}__recon_mni/\${bname} \
        --keep_invalid
    done
    """
}

bundles_for_longitudinal
    .groupTuple(by: 0)
    .map{[it[0], it[2]]}
    .set{atlas_for_longitudinal}

tractograms_for_recognition_second
    .join(anat_for_reference_bundles_second, by: [0,1])
    .join(transformation_for_recognition_second, by: [0,1])
    .combine(atlas_config_second)
    .combine(atlas_for_longitudinal, by: 0)
    .set{tractograms_and_transformation_individual}
process Recognize_Bundles_Individual {
    publishDir = {"./results_rbx/$sid/$sess/$task.process/"}
    memory { params.single_dataset_size_GB.GB * params.rbx_processes }
    cpus params.rbx_processes

    input:
    set sid, sess, file(tractograms), file(reference), file(transfo),
        file(config), file(directories) from tractograms_and_transformation_individual
    output:
    set sid, sess, "*.trk" into bundles_for_cleaning
    script:
    """
    mkdir tmp/
    mkdir atlas/
    mv ${directories} atlas/
    scil_recognize_multi_bundles.py $tractograms ${config} atlas/ ${transfo} \
        --inverse --out_dir tmp/ \
        --log_level DEBUG --minimal_vote_ratio $params.minimal_vote_ratio \
        --seed $params.seed --processes $params.rbx_processes
    mv tmp/* ./
    """
}

bundles_for_cleaning
    .transpose()
    .set{all_bundles_for_cleaning}
process Clean_Bundles {
    publishDir = {"./results_rbx/$sid/$sess/$task.process/"}

    input:
    set sid, sess, file(bundle) from all_bundles_for_cleaning

    output:
    set sid, sess, val(bname), "${sid}__*_cleaned.trk" optional true into bundle_for_density
    script:
    bname = bundle.name.take(bundle.name.lastIndexOf('.'))
    """
    scil_outlier_rejection.py ${bundle} "${sid}__${bname}_cleaned.trk" \
        --alpha $params.outlier_alpha
    """
}

bundle_for_density
    .combine(transformation_for_average, by: [0, 1])
    .combine(atlas_anat_for_average)
    .set{all_bundles_transfo_for_average}
process Compute_Density_Bundles {
    publishDir = {"./results_rbx/$sid/$sess/$task.process/"}

    input:
    set sid, sess, val(bname), file(bundle), file(transfo), file(atlas) from all_bundles_transfo_for_average
    
    output:
    set bname, "${sid}-${sess}__${bname}_density.nii.gz",
        "${sid}-${sess}__${bname}_binary.nii.gz" into bundle_for_average
    when:
    params.run_average_bundles
    script:
    """
    scil_apply_transform_to_tractogram.py $bundle $atlas $transfo tmp.trk \
        --remove_invalid
    scil_compute_streamlines_density_map.py tmp.trk "${sid}-${sess}__${bname}_density.nii.gz"
    scil_image_math.py lower_threshold "${sid}-${sess}__${bname}_density.nii.gz" 1 \
        "${sid}-${sess}__${bname}_binary.nii.gz"
    """
}


bundle_for_average
    .groupTuple(by: 0)
    .set{all_bundle_for_average}
process Average_Bundles {
    publishDir = params.Average_Bundles_Publish_Dir

    input:
    set val(bname), file(bundles_den), file(bundles_bin) from all_bundle_for_average

    output:
    file "${bname}_density.nii.gz"
    file "${bname}_binary.nii.gz"
    script:
    """
    scil_image_math.py addition ${bundles_den} 0 ${bname}_density.nii.gz
    scil_image_math.py addition ${bundles_bin} 0 ${bname}_binary.nii.gz
    """
}
