source(paste0(config$cobase_dir, "Postprocessing/Utilities/uvpp_util.R"))
source(paste0(config$cobase_dir, "Postprocessing/Utilities/sim_util.R"))
source(paste0(config$cobase_dir, "Postprocessing/Utilities/ensfc_util.R"))




preprocess_all <- function(file_name, observation_columns, ensemble_regex){
    # Load the csv data
    data <- load_data(file_name)

    # Apply UVPP and save
    uvpp <- compute_uvpp(data = data, 
                         group_name = file_name, 
                         observation_columns = observation_columns, 
                         ensemble_regex = ensemble_regex)

    # Transform the data and save
    transformed_data <- transform_data(data, file_name, observation_columns, ensemble_regex)

    # Compute similarity matrix
    compute_sim_matrix(data, file_name)
}

# Mock data
preprocess_all(
    file_name               = "Mock_data",
    observation_columns     = c("obs"),
    ensemble_regex          = c("^M_")
)
