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
Channel
    .fromPath("$root/**/*tracking*.trk",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .set{tractogram_for_recognition}

Channel
    .fromPath("$root/**/*fa.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{anat_for_registration;anat_for_reference_centroids;anat_for_reference_bundles}

if (!(params.atlas_anat) || !(params.atlas_config) || !(params.atlas_directory)) {
    error "You must specify all 3 atlas related input. --atlas_anat, " +
    "--atlas_config and --atlas_directory all are mandatory."
}

Channel.fromPath("$params.atlas_anat")
    .into{atlas_anat;atlas_anat_for_average}
atlas_config = Channel.fromPath("$params.atlas_config")
atlas_directory = Channel.fromPath("$params.atlas_directory")

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


tractogram_for_recognition
    .join(anat_for_reference_bundles)
    .join(transformation_for_recognition)
    .combine(atlas_config)
    .combine(atlas_directory)
    .set{tractogram_and_transformation}
process Recognize_Bundles {
    cpus params.rbx_processes
    memory params.rbx_memory_limit

    input:
    set sid, file(tractogram), file(refenrence), file(transfo), file(config), file(directory) from tractogram_and_transformation
    output:
    set sid, "*.trk" into bundles_for_cleaning
    file "results.json"
    file "logfile.txt"
    script:
    """
    scil_remove_invalid_streamlines.py ${tractogram} tractogram_ic.trk --reference ${refenrence} --remove_single_point --remove_overlapping_points
    mkdir tmp/
    scil_recognize_multi_bundles.py tractogram_ic.trk ${config} ${directory}/*/ ${transfo} --inverse --out_dir tmp/ \
        --log_level DEBUG --multi_parameters $params.multi_parameters --minimal_vote_ratio $params.minimal_vote_ratio \
        --tractogram_clustering_thr $params.wb_clustering_thr --seeds $params.seeds --processes $params.rbx_processes
    rm tractogram_ic.trk
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
    publishDir = "./results_rbx/Average_Bundles/"
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
