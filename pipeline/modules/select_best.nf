process SELECT_BEST {
    tag "best_parameters"
    label 'low'
    container "docker://mdagnelie/parameter_estimation:v0.1"
    publishDir "${params.outdir}/plots", mode: 'copy'
    
    input:
    path(all_fits)
    path(dataset)
    
    output:
    path "best_phi.json"
    path "objective_plot.png"
    path "best_fit.png"
    
    script:
    """
    select_best.py ${all_fits} ${dataset}
    """
}
