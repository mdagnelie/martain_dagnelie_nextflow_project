include {CSVTK_CONCAT} from './modules/nf-core/csvtk/concat/main.nf'
include {GENERATE_GUESS} from './modules/generate_guess.nf'
include {SIMULATION_FIT} from "./modules/simulation_fit.nf"
include {SELECT_BEST} from "./modules/select_best.nf"

workflow {
    //combine the time series
    def data_files_ch = Channel.fromPath(params.time_series)
        .collect()
        .map { time_serie -> tuple([id:"all_timeseries"], time_serie)}

    dataset_ch= CSVTK_CONCAT(data_files_ch,"csv","csv").csv 
        .map{ meta,file -> file }

    guesses_ch = GENERATE_GUESS(params.n_samples).flatten()
    
    combined_results = SIMULATION_FIT(guesses_ch, dataset_ch)
        .collect()                    
        .flatten()                    
        .collectFile(name: 'all_fits.csv', storeDir: "${params.outdir}/samples_fit") 
    
    SELECT_BEST(combined_results, dataset_ch.first())

}


