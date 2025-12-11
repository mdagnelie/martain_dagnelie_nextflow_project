process SIMULATION_FIT {
    tag 'likelihood'
    label 'high'
    container "docker://mdagnelie/parameter_estimation:v0.1"
    
    input:
    path(phi_json)
    path(dataset)
    
    output:
    path 'fit_result.csv'
    
    script:
    """
    simulation_fit.py ${phi_json} ${dataset}
    """
}