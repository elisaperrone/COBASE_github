source(paste0(config$cobase_dir, "/Postprocessing/Utilities/mvpp.R"))

# Mock data
postprocess_all(
    file_name               = "Mock_data",
    observation_columns     = c("obs"),
    ensemble_regex          = c("^M_"),
    output_dim_standard     = 10
)
