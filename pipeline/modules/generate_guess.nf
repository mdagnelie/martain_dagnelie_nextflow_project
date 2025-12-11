process GENERATE_GUESS{
    tag "generate_guess"
    label 'low'
    container "docker://mdagnelie/parameter_estimation:v0.1"

    input:
    val(n_samples)

    output:
    path "*.json"
 
    script:
    """
    generate_guess.py ${n_samples}
    """
}