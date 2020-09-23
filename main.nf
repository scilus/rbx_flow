#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["atlas_config":"$params.atlas_config",
                "atlas_directory":"$params.atlas_directory",
                "atlas_centroids":"$params.atlas_centroids",
                "multi_parameters":"$params.multi_parameters",
                "minimal_vote_ratio":"$params.minimal_vote_ratio",
                "wb_clustering_thr":"$params.wb_clustering_thr",
                "seeds":"$params.seeds",
                "outlier_alpha":"$params.outlier_alpha",
                "register_processes":"$params.register_processes",
                "rbx_processes":"$params.rbx_processes"]

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

log.info "Input: $params.root"
root = file(params.root)
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
    .into{anat_for_registration;anat_for_reference_c;anat_for_reference_b}

atlas_anat = Channel.fromPath("$params.atlas_anat")
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
    input:
    set sid, file(native_anat), file(atlas) from anats_for_registration

    output:
    set sid, "${sid}__output0GenericAffine.mat" into transformation_for_recognition, transformation_for_centroids
    file "${sid}__outputWarped.nii.gz"
    script:
    """
    antsRegistrationSyNQuick.sh -d 3 -f ${native_anat} -m ${atlas} -n ${params.register_processes} -o ${sid}__output -t a
    """ 
}


anat_for_reference_c
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
    .join(anat_for_reference_b)
    .join(transformation_for_recognition)
    .combine(atlas_config)
    .combine(atlas_directory)
    .set{tractogram_and_transformation}
process Recognize_Bundles {
    cpus params.rbx_processes
    input:
    set sid, file(tractogram), file(refenrence), file(transfo), file(config), file(directory) from tractogram_and_transformation
    output:
    set sid, "*.trk" into bundles_for_cleaning
    file "results.json"
    file "logfile.txt"
    script:
    """
    scil_remove_invalid_streamlines.py ${tractogram} tractogram_ic.trk --reference ${refenrence}
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
    file "${sid}__*_cleaned.trk"
    script:
    bname = bundle.name.take(bundle.name.lastIndexOf('.'))
    """
    scil_outlier_rejection.py ${bundle} "${sid}__${bname}_cleaned.trk" --alpha $params.outlier_alpha
    """ 
}
