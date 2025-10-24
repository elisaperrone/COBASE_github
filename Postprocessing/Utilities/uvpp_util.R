compute_uvpp <- function(data, group_name, observation_columns, ensemble_regex, ai_regex = NULL) {
    # Data structure to store all data
    uvpp <- list()

    # UVPP is done for each station and variable
    stations <- unique(data$station)

    cat("In Postprocessing/uvpp_util.R that is called in preprocess.R: Computing UVPP output...\n")

    # Apply UVPP per (station, obs_column) group
    for (var_idx in 1:length(observation_columns)) {
        for (station_nr in stations) {
            # Observation information
            obs_column <- observation_columns[var_idx]
            ens_regex <- ensemble_regex[var_idx]

            # Diagnostic information
            cat(paste0(
                format(Sys.time(), "%H:%M:%S"),
                " | In Postprocessing/uvpp_util.R that is called in univariate.R: Computing UVPP for (station, obs) = (", station_nr, ", ", obs_column, ")\n"
            ))

            # Subset on station
            uvpp_data <- subset(data, station == station_nr)

            if (!is.null(ai_regex)) {
                ai_reg <- ai_regex[var_idx]

                # EMOS computation
                emos.result <- emos_T2M_mean_singleForecast_IFS_and_AI(
                    uvpp_data,
                    config$UVPP_Training_Window,
                    ens_regex,
                    ai_reg,
                    obs_column
                )
            } else {
                # EMOS computation
                emos.result <- emos_T2M_mean_singleForecast(
                    uvpp_data,
                    config$UVPP_Training_Window,
                    ens_regex,
                    obs_column
                )
            }



            # Save the results
            emos.result$station <- station_nr
            emos.result$obs_type <- obs_column
            uvpp <- rbind(uvpp, emos.result)
        }
    }

    cat("\n")

    savename <- paste0(config$cobase_dir, config$results_folder, "/UVPP/uvpp_", group_name, ".Rdata")
    if (!dir.exists(dirname(savename))) {
        dir.create(dirname(savename), recursive = TRUE)
    }
    cat(paste0(
      format(Sys.time(), "%H:%M:%S"),
      " | In Postprocessing/uvpp_util.R that is called in univariate.R: saving uvpp to", savename, "\n"
    ))
    
    save(uvpp, file = savename)

    return(uvpp)
}


# -------------------------------------------------------------------
# - NAME:        ECC_T2M_Emos_subfunctions.R
# - AUTHOR:      Perrone Elisa
# - DATE:        2021-02-18
# -------------------------------------------------------------------
# - DESCRIPTION: EMOS -- optimization based on CRPS
# -------------------------------------------------------------------
# - CONTENT:    -> emos_T2M_mean_singleForecast
#                  (tailored for ALADIN-LAEF data)
# -------------------------------------------------------------------


#####################################################################
#  Function that applies EMOS to a dataset data.lt
#  leadtime and stations are fixed
#  (adjusted from Moritz's code)
#####################################################################

emos_T2M_mean_singleForecast <- function(
    data.lt,
    trainingDays,
    ens_regex,
    obs_column) {
    # Group the ensemble members
    ens_data <- data.lt[, grep(ens_regex, names(data.lt))]

    # Ensemble statistics
    data.lt$ensmean <- apply(ens_data, 1, mean)
    data.lt$enssd <- apply(ens_data, 1, sd)
    data.lt$ensvar <- apply(ens_data, 1, var)
    data.lt$ensmedian <- apply(ens_data, 1, median)

    nd <- length(unique(data.lt$td))
    td.emos <- trainingDays # training days
    npred <- 1

    # model estimation
    slopeCoef <- array(dim = c((nd - td.emos), (npred + 1)))
    varCoef <- array(dim = c((nd - td.emos), 2))
    valid_df <- data.frame(matrix(nrow = (nd - td.emos), ncol = 6))

    # loop over days
    for (i.d in 1:(nd - td.emos)) {
        train.index <- data.lt$ytime %in% data.lt$ytime[i.d:(i.d + td.emos - 1)]
        valid.index <- data.lt$ytime %in% data.lt$ytime[(i.d + td.emos)]

        # fit linear regression for initial values
        lm.fit <- lm(as.formula(paste(obs_column, "~ ensmean")), data = data.lt[train.index, ])

        # (3b) use initial values according to R Package of Gneiting et al. 2005
        a.init <- coef(lm.fit)[1]
        B.init <- coef(lm.fit)[-1]
        B.init <- sqrt(abs(B.init))
        c.init <- 5
        c.init <- sqrt(c.init)
        d.init <- 1
        d.init <- sqrt(d.init)

        par_init <- c(a.init, B.init, c.init, d.init)

        t.obs <- data.lt[train.index, obs_column]
        t.ens <- cbind(1, as.matrix(data.lt$ensmean[train.index]))
        t.var <- data.lt$ensvar[train.index]


        # (3c) fitting of parameters using CRPS minimization
        fits <- optim(par_init, crps.normal,
            t.obs = t.obs,
            t.ens = t.ens,
            t.var = t.var,
            K = dim(t.ens)[2],
            method = "BFGS",
            control = list(maxit = as.integer(1e7))
        )

        #  (3d) transform parameters to observable scale and store them in data.frame
        fits$par[-1] <- fits$par[-1]^2
        slopeCoef[i.d, ] <- fits$par[1:(npred + 1)]
        varCoef[i.d, ] <- fits$par[(npred + 2):length(fits$par)]

        #  (4) VALIDATION
        #  (4a) calculate mean, variance, standard deviation
        ens_mu <- slopeCoef[i.d, ] %*% c(1, data.lt$ensmean[valid.index])
        ens_var <- varCoef[i.d, ] %*% c(1, data.lt$ensvar[valid.index])
        ens_sd <- sqrt(ens_var)

        #  (4b) calculate crps
        crps.emos <- crps.norm(data.lt[valid.index, obs_column], ens_mu, ens_sd)
        crps.raw_ens <- crps.ensemble(
            as.numeric(data.lt[valid.index, obs_column]),
            as.matrix(data.lt[valid.index, grep(ens_regex, names(data.lt))])
        )

        #  (4c)  store scores in data.frame
        valid_df[i.d, ] <- cbind(data.lt$td[valid.index], data.lt[valid.index, obs_column], ens_mu, ens_sd, crps.emos, crps.raw_ens)
        names(valid_df) <- c("td", "obs", "ens_mu", "ens_sd", "crps_emos", "crps.raw_ens")
    } # end of loop over days

    return(valid_df)
}


# -------------------------------------------------------------------
# - NAME:        emos_T2M_mean_singleForecast_subfunctions.R
# - AUTHOR:      Moritz Lang
# - DATE:        2016-07-20
# -------------------------------------------------------------------
# - DESCRIPTION: contains subfunctions for emos_T2M_mean_singleForecast.R
# -------------------------------------------------------------------
# - CONTENT:  - crps.normal()
#             - crps.ensemble()
#             - crps.norm()
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# CRPS calculation as in R-package
# (adapted from Gneiting et al. 2005)
# -------------------------------------------------------------------
crps.normal <- function(param, t.obs, t.ens, t.var, K) {
    param[-1] <- param[-1]^2
    MEAN <- t.ens %*% param[1:K]
    VAR <- param[K + 1] + param[K + 2] * t.var
    SIGMA <- sqrt(VAR)
    obs.stdz <- (t.obs - MEAN) / SIGMA
    crps.obs <- SIGMA * (obs.stdz * (2 * .Internal(pnorm(obs.stdz, 0, 1, TRUE, FALSE)) - 1) + 2 * .Internal(dnorm(obs.stdz, 0, 1, FALSE)) - 1 / sqrt(pi))

    return(.Primitive("sum")(crps.obs)) # mean to sum changed!!
}
library("compiler")
crps.normal <- cmpfun(crps.normal)

# -------------------------------------------------------------------
# Continuous Rank Probability Score: on discret number of forecasts, ensemble
# (adapted from Reto Stauffer 2016)
# -------------------------------------------------------------------
crps.ensemble <- function(obs, fcst, plot = FALSE) {
    # Sorting ensemble
    fcst <- sort(fcst)
    # Ensemble empirical distribution
    fun.ecdf <- ecdf(fcst)
    # Differences and corresponding cdf's
    diff <- diff(sort(c(fcst, obs)))
    cdfs <- fun.ecdf(sort(c(fcst, obs)))[1:length(fcst)]
    # Compute crps
    crps <- sum((cdfs - (fcst > obs))^2 * diff)
    # Demo plot
    if (plot) {
        plot(fcst, fun.ecdf(fcst), type = "s", main = sprintf(" RPS score: %.3f", crps))
        abline(v = obs, col = 2)
    }
    crps
}

# -------------------------------------------------------------------
# CRPS for normal normal distribution
# (adapted from Gneiting et al. 2005)
# -------------------------------------------------------------------
crps.norm <- function(obs, mu, sd) {
    sd * ((obs - mu) / sd * (2 * pnorm((obs - mu) / sd) - 1) + 2 * dnorm((obs - mu) / sd) - 1 / sqrt(pi))
}

crps.normal_split <- function(param, t.obs, t.ens, t.var, K) {
    param[-1] <- param[-1]^2
    MEAN <- t.ens %*% param[1:K]
    VAR <- param[K + 1] + param[K + 2] * t.var[, 1] + param[K + 3] * t.var[, 2]
    SIGMA <- sqrt(VAR)
    obs.stdz <- (t.obs - MEAN) / SIGMA
    crps.obs <- SIGMA * (obs.stdz * (2 * pnorm(obs.stdz, 0, 1) - 1) + 2 * dnorm(obs.stdz, 0, 1) - 1 / sqrt(pi))
    sum(crps.obs)
}
crps.normal_split <- cmpfun(crps.normal_split)
