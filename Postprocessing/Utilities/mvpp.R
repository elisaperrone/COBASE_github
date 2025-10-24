source(paste0(config$cobase_dir, "/Postprocessing/Utilities/load_data_util.R"))
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/ensfc_util.R"))
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/scores_util.R"))
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/mvpp_methods.R"))
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/DM_util.R"))

postprocess_all <- function(file_name, observation_columns, ensemble_regex, output_dim_standard)
{

    # Fixed settings
    trainingDays            <- config$UVPP_Training_Window + config$MVPP_Training_Window
    trainingWindow          <- config$MVPP_Training_Window

    # set random seed
    set.seed(1)

    # Stations and days from the data
    data            <- load_data(file_name)
    stations        <- unique(data$station)
    days            <- sort(unique(data$td))
    nout            <- length(days) - trainingDays
    d               <- length(stations) * length(observation_columns)
    get_dimension   <- load_dimension_transform(stations, observation_columns)

    # Ensemble members
    m <- sum(grepl(ensemble_regex[1], names(data))) # Each observation should have the same number of ens members

    # Environment to store scores in
    score.env <- new.env()

    # Save settings
    score.env$file_name             <- file_name
    score.env$d                     <- d
    score.env$stations              <- as.character(stations) # keep the station names so I can find them later
    score.env$nout                  <- nout
    score.env$trainingDays          <- trainingDays
    score.env$trainingWindow        <- trainingWindow
    score.env$observation_columns   <- observation_columns
    score.env$ensemble_regex        <- ensemble_regex

    # Load transformed data
    transformed_data                <- load_transformed_data(file_name)
    score.env$obs                   <- transformed_data$obs
    score.env$obs_init              <- transformed_data$obs
    score.env$ensfc                 <- transformed_data$ensfc
    score.env$ensfc_init            <- transformed_data$ensfc_init
    score.env$mvpp_list$ens         <- transformed_data$ensfc

    # Load precomputed datastructures
    sim_matrix                      <- load_sim_matrix(file_name)
    score.env$sim_matrix            <- sim_matrix

    # Add the ens scores
    add_scores(score.env = score.env,
               res = list(mvppout = transformed_data$ensfc), 
               method = "ens",
               start_time = 0)

    # General function to apply a MVPP method
    mvpp_adjusted <- function(...) {
        if ("output_dim" %in% names(list(...))) {
            return(mvpp(...,
                transformed_data  = transformed_data,
                score.env         = score.env
            ))
        }

        if ("trainingWindow" %in% names(list(...))) {
            return(mvpp(...,
                transformed_data  = transformed_data,
                score.env         = score.env,
                output_dim        = output_dim_standard,
                addTrainingDays   = TRUE
            ))
        }

        return(mvpp(...,
            transformed_data  = transformed_data,
            score.env         = score.env,
            trainingWindow    = trainingWindow,
            output_dim        = output_dim_standard
        ))
    }



    ##################
    ## EMOS for ECC ##
    ##################
    print("In Prosprocessing/Utilities/mvpp.R: doing the post-processing...")

    emos.q <- mvpp_adjusted(
        method            = "EMOS",
        variant           = "Q",
        output_dim        = m,
        saveScores        = FALSE
    )

    emos.r <- mvpp_adjusted(
        method            = "EMOS",
        variant           = "R",
        output_dim        = m,
        saveScores        = FALSE
    )

    #########
    ## ECC ##
    #########

    ecc.q <- mvpp_adjusted(
        method            = "ECC",
        variant           = "Q", # For saved name
        EMOS_sample       = emos.q$mvppout
    )

    ecc.r <- mvpp_adjusted(
        method            = "ECC",
        variant           = "R", # For saved name
        EMOS_sample       = emos.r$mvppout
    )

    ##########
    ## EMOS ##
    ##########

    emos.q <- mvpp_adjusted(
        method            = "EMOS",
        variant           = "Q"
    )

    emos.r <- mvpp_adjusted(
        method            = "EMOS",
        variant           = "R"
    )

    #########
    ## SSH ##
    #########
    # Regular Schaake Shuffle
    ssh.i <- mvpp_adjusted(
        method            = "SSh-I14",
        EMOS_sample       = emos.q$mvppout,
        variant           = "Q"
    )

    # Sim Schaake 
    ssh.sim <- mvpp_adjusted(
      method            = "SimSchaake", 
      EMOS_sample       = emos.q$mvppout,
      sim_matrix        = sim_matrix,
      variant           = "Q"
    )
    
    #Sim Schaake with all past observations - version R
    ssh.sim <- mvpp_adjusted(
      method            = "SimSchaake", 
      EMOS_sample       = emos.r$mvppout,
      sim_matrix        = sim_matrix,
      variant           = "R"
    )

    #########
    ## GCA ##
    #########

    gca.cop <- mvpp_adjusted(
        method = "CopGCA"
    )

    # EMOS-Q Shuffling
    gca.cop.sh <- mvpp_adjusted(
        method            = "CopGCA",
        EMOS_sample       = emos.q$mvppout,
        shuffle           = TRUE,
        MVPP_sample       = gca.cop$mvppout,
        variant           = "Q"
    )

    #########################
    ## Archimedean Copulas ##
    #########################

    clayton <- mvpp_adjusted(
        method = "Clayton"
    )

    gumbel <- mvpp_adjusted(
        method = "Gumbel"
    )

    frank <- mvpp_adjusted(
        method = "Frank"
    )

    # EMOS-Q Shuffling
    clayton.sh <- mvpp_adjusted(
        method            = "Clayton",
        EMOS_sample       = emos.q$mvppout,
        shuffle           = TRUE,
        MVPP_sample       = clayton$mvppout,
        variant           = "Q"
    )

    gumbel.sh <- mvpp_adjusted(
        method            = "Gumbel",
        EMOS_sample       = emos.q$mvppout,
        shuffle           = TRUE,
        MVPP_sample       = gumbel$mvppout,
        variant           = "Q"
    )

    frank.sh <- mvpp_adjusted(
        method            = "Frank",
        EMOS_sample       = emos.q$mvppout,
        shuffle           = TRUE,
        MVPP_sample       = frank$mvppout,
        variant           = "Q"
    )

    #####################
    ## Save the scores ##
    #####################
    # putting it into the folder results / version / scores 
    score_dir <- paste0(config$cobase_dir, config$results_folder, "/Scores/") 

    if (!dir.exists(score_dir)) {
        dir.create(score_dir, recursive = TRUE)
    }

    savename <- paste0(file_name, "_mout_", output_dim_standard)

    save(
        list = ls(score.env),
        file = paste0(score_dir, "score_env_", savename, ".RData"),
        envir = score.env
    )

    ###########################
    ## Compute DM statistics ##
    ###########################
    #config$all_benchmarks 
    print("In Postprocessing/Utilities/mvpp.R: computing DM scores...")
    benchmarks = c(config$dmcrps_bm, config$dmmvsimssh_bm, config$dmesvs_mvpp_bm, config$dmesvs_parcop_bm)
    compute_DM_scores(file_name = savename, 
                      benchmarks = benchmarks, 
                      parallelization = FALSE)

}
