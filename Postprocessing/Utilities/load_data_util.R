library(dplyr)

load_data <- function(file_name = "ifs_and_kis_dataset") {
    # Load KNMI observations from KIS (Klimatologisch Informatie Systeem)
    all_data <- read.csv(paste0(config$cobase_dir, config$data_folder, "/Datasets/", file_name, ".csv"))

    # Transformations to help R understand the data
    all_data$date <- as.Date(all_data$validTime)
    all_data$station <- factor(all_data$station)

    # Remove unused columns
    all_data$X <- NULL
    all_data$validTime <- NULL
    all_data$runtime <- NULL

    # Add date indices
    all_data <- all_data %>%
        arrange(date) %>%
        mutate(td = dense_rank(date))

    # For compatibility with other functions
    all_data$ytime <- all_data$date

    return(all_data)
}

load_uvpp <- function(file_name) {
    # Loads a variable called uvpp
    load(paste0(config$cobase_dir, config$results_folder, "/UVPP/uvpp_", file_name, ".Rdata"))

    return(uvpp)
}

load_sim_matrix <- function(file_name) {
    # Loads a variable called uvpp
    load(paste0(config$cobase_dir, config$results_folder, "/SimilarityMatrix/simMatrix_", file_name, ".Rdata"))

    return(sim_matrix)
}

load_dimension_transform <- function(stations, observation_columns) {
    get_dimension <- function(stat, obs) {
        stat_idx <- match(stat, stations)
        obs_idx <- match(obs, observation_columns)

        return(stat_idx + length(stations) * (obs_idx - 1))
    }

    return(get_dimension)
}

load_transformed_data <- function(file_name) {
    # Loads a variable called transformed_data
    load(paste0(config$cobase_dir, config$data_folder, "/ENS/transformed_data_", file_name, ".Rdata"))

    return(transformed_data)
}

load_score_env <- function(file_name, envir) {
    fName <- paste0(config$cobase_dir, config$results_folder, "/Scores/score_env_", file_name, ".RData")
    cat("In Postprocessing/Utilities/load_data_util.R: fName == ", fName,"\n")
    if (file.exists(fName)) {
        # Loads many score variables
        load(fName, envir = envir)
    } else 
    {
        cat("In Postprocessing/Utilities/load_data_util.R: The file", fName, "does not exist.\n")
    }
}

load_scores <- function(res) {
    input_scores <- c(
        "crps_list",
        "es_list",
        "vs0_list",
        "vs0w_list",
        "vs1_list",
        "vs1w_list"
    )
    d <- res$d

    # this is one place where the names of crps_1 comes from
    #for (dd in 1:d)
    #{
    #    input_scores <- c(input_scores, paste0("crps_", dd))
    #}
    input_scores = c(input_scores, paste0("crps_", as.vector(outer(res$observation_columns, res$stations, paste, sep = "-"))))

    return(input_scores)
}

load_dm_statistics <- function(file_name = "all_data_kiri", comparison = FALSE) {
    # Loads a variable
    if (comparison) {
        load(paste0(config$cobase_dir, config$results_folder, "/TestStatistic/dm_statistics_comparison_", file_name, ".Rdata"))
    } else {
        load(paste0(config$cobase_dir, config$results_folder, "/TestStatistic/dm_statistics_", file_name, ".Rdata"))
    }

    return(all_dfmc)
}

pretty_score_name <- function(score_name) {
    if (grepl("^crps_\\d+$", score_name)) {
        return(sub("^crps_(\\d+)$", "CRPS-\\1", score_name))
    }

    pretty_name <- switch(score_name,
        "crps_list"         = "CRPS",
        "es_list"           = "ES",
        "vs0_list"          = "VS-0.5",
        "vs0w_list"         = "VS-w-0.5",
        "vs1_list"          = "VS-1",
        "vs1w_list"         = "VS-w-1",
        "???"
    )

    return(pretty_name)
}

pretty_model_name <- function(model_name) {
    if (model_name == "SimSchaake-CopGCAsh") return ("SimShuffled GCA")

    if (grepl("-66sh$", model_name)) {
        return(paste0(pretty_model_name(sub("-66sh$", "", model_name)), "-66-sh"))
    } else if (grepl("sh$", model_name)) {
        return(paste0(pretty_model_name(sub("sh$", "", model_name)), "-sh"))
    } else if (grepl("^Sim", model_name)) {
        return(paste0("Sim-", pretty_model_name(sub("^Sim", "", model_name))))
    }

    pretty_name <- switch(model_name,
        "ens"           = "Raw Ensemble",
        "CopGCA"        = "GCA",
        "SimSchaake-H"  = "SimSchaake",
        model_name # Not found
    )

    return(pretty_name)
}
