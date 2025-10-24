source(paste0(config$cobase_dir, "/Creating Plots/Utilities/dm_plots.R")) # defines the function: load_and_prepare_data
source(paste0(config$cobase_dir, "/Postprocessing/Utilities/load_data_util.R")) # defines the function: load_score_env

library(xtable)


########################################
# load functions
########################################
change_names <- function(model_name) {
  pretty_name <- switch(model_name,
                        "ens"           = "Raw Ensemble",
                        "CopGCA"        = "GCA",
                        "SimSchaake-Q" = "SimSSh",
                        "SimSchaake-R" = "SimSSh-R",
                        "ECC-Q" = "ECC",
                        "SSh-I14-Q" = "SSh",
                        "COBASE-CopGCA-Q" = "COBASE-GCA",
                        "COBASE-Clayton-Q" = "COBASE-Clayton",
                        "COBASE-Frank-Q" = "COBASE-Frank",
                        "COBASE-Gumbel-Q" = "COBASE-Gumbel",
                        
                        model_name # Not found
  )
  
  return(pretty_name)
}

########################################
# define common info and vars
########################################
# https://www.knmi.nl/nederland-nu/klimatologie/daggegevens
# https://dataset.api.hub.geosphere.at/app/frontend/station/historical/synop-v1-1h
# https://acinn-data.uibk.ac.at/pages/tawes-uibk.html

station_info = data.frame(st_num = c("1", "2", "3"),
                          st_name = c("St1", "St2", "St3"),
                          group = c("Mock_data","Mock_data", "Mock_data"),
                          lat = c(NA, NA, NA),
                          lon = c(NA, NA, NA))


all_groups_info = data.frame(group_names = "Mock_data",
                             group_names_pretty = "Mock Data",
                             output_dim_standard = 10)

# alpha for plotting
alpha <- 0.25

# 
########################################
########## dmcrps_t2m: (Figure 1)
########################################
# CRPS EMOS-R v EMOS-Q
print(paste0("working on: dmcrps_t2m"))
dmcrps_benchmark_name  = config$dmcrps_bm
dmcrps_model_names     = c(config$dmcrps_models)

# loop over groups:
dmcrps_scores = lapply(seq_along(all_groups_info$group_names), function(gn){
  print(paste0("gn == ", gn))
  # get group name:
  file_name  <- paste0(all_groups_info$group_names[gn], "_mout_", all_groups_info$output_dim_standard[gn])


  # load scores:
  res <- new.env()
  load_score_env(file_name, res)

  # get the names of the scores that we want to keep:
  crps_scores <- paste0("crps_", as.vector(outer(res$observation_columns, res$stations, paste, sep = "-")))

  # load the actual data with the bootstrap samples of all the scores:
  data_env <- load_and_prepare_data(file_name, dmcrps_benchmark_name)
  df <- data_env$df
  # filter by the scores that we want to keep:
  dfplot <- filter_models(df, dmcrps_model_names, crps_scores)$df %>% subset(score %in% crps_scores)
  # reformat for plotting:
  dfplot <- dfplot %>%
    # pivot
    pivot_longer(cols = starts_with("bootstrap_"), names_to = "bootstrap", values_to = "value") %>%
    # assign a station number from score name, and rename label the score to crps
    mutate(station = gsub("crps_", "", score), score = "crps") %>%
    # separate the variable and the station number
    separate(station, into = c("var", "station"), sep = "-") %>%
    #recode the variable for plotting
    mutate(var = ifelse(var %in% c("obs", "T_DRYB_10"), "T2m", "DPT"))

  dfplot$group = all_groups_info$group_names[gn]
  return(dfplot)
  }) %>% bind_rows()


# join with the pretty station info
dmcrps_scores <- dmcrps_scores %>%
  left_join(station_info, by = c("station" = "st_num"))

# put the stations in order
dmcrps_scores$st_name = factor(dmcrps_scores$st_name, levels = station_info$st_name)


pdf(file = paste0(config$cobase_dir, config$results_folder, "Figures/flos_dmcrps_t2m.pdf"), width = 12, height = 6)
print(
ggplot(dmcrps_scores %>% filter(var == "T2m"), aes(st_name, value)) +
  facet_wrap(~var) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = qnorm(alpha), ymax = qnorm(1 - alpha)),
            fill = "gray75", color = "gray75", alpha = alpha) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray25") +
  theme_bw() +
  xlab("Station") + ylab("DM test statistic") +
  ggtitle(label = "CRPS: EMOS v EMOS-R") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(size=16),
        plot.margin = margin(t = 10, r = 33, b = 10, l = 45))
)
dev.off()

########################################
########## dmmvsimssh: (Figure 2)
########################################
#
print(paste0("working on: dmmvsimssh"))

dmmvsimssh_model_names     = c(config$dmmvsimssh_models)
dmmvsimssh_benchmark_name  = config$dmmvsimssh_bm
dmmvsimssh_scores_list = c("es_list", "vs1_list")


# loop over groups
dmmvsimssh_scores = lapply(seq_along(all_groups_info$group_names), function(gn){
  print(paste0("gn == ", gn))
  # get group name:
  file_name  <- paste0(all_groups_info$group_names[gn], "_mout_", all_groups_info$output_dim_standard[gn])

  # loop over the energy and variogram scores
  dmmvsimssh_scores_tmp = list()
  for(sc in seq_along(dmmvsimssh_scores_list)){
    keep_scores = dmmvsimssh_scores_list[sc] #"es_list"

    # load the actual data and get the data.frame, and filter:
    data_env <- load_and_prepare_data(file_name, dmmvsimssh_benchmark_name)
    df <- data_env$df
    dfplot <- filter_models(df, dmmvsimssh_model_names, keep_scores)$df %>% subset(score %in% keep_scores)
    # rearrange for plotting
    dfplot <- dfplot %>%
      pivot_longer(cols = starts_with("bootstrap_"), names_to = "bootstrap", values_to = "value") %>%
      # recode some things for the plotting
      mutate(group_name = file_name,
             group_names_pretty = all_groups_info$group_names_pretty[gn],
             score = ifelse(dmmvsimssh_scores_list[sc] == "es_list", "Energy Score", "Variogram Score"))
    # get prettier names for the models
    new_model_tmp = unlist(lapply(as.character(dfplot$model), change_names))
    dfplot$model = new_model_tmp
    # get prettier names for the bm
    new_bm_tmp = unlist(lapply(as.character(dfplot$benchmark), change_names))
    dfplot$benchmark = new_bm_tmp

    # save for this score
    dmmvsimssh_scores_tmp[[sc]] = dfplot

  }
    return(bind_rows(dmmvsimssh_scores_tmp) )
}) %>% bind_rows()


# get prettier group names for the plot
dmmvsimssh_scores$group_names_pretty = factor(gsub("_mout_[0-9][0-9]", "", dmmvsimssh_scores$group_names_pretty), levels = all_groups_info$group_names_pretty)


pdf(file = paste0(config$cobase_dir, config$results_folder, "Figures/flos_dmesvs_simssh.pdf"), width = 12, height = 6)
print(
ggplot(dmmvsimssh_scores, aes(group_names_pretty, value)) +
  facet_wrap(~score, ncol = 2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = qnorm(alpha), ymax = qnorm(1 - alpha)),
            fill = "gray75", color = "gray75", alpha = alpha) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray25") +
  theme_bw() +
  xlab("Data group") + ylab("DM test statistic") +
  ggtitle(label = paste0("SimSSh v ", "SimSSh-R")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(size=16),
        plot.margin = margin(t = 10, r = 33, b = 10, l = 45))
)
dev.off()
  
  
  
########################################
########## dmesvs_mvpp (Figure 3):
########################################
# 
print(paste0("working on: dmesvs_mvpp"))

dmesvs_mvpp_model_names     = c(config$dmesvs_mvpp_models)
dmesvs_mvpp_benchmark_name  = config$dmesvs_mvpp_bm[2] 
dmesvs_mvpp_scores_list = c("es_list", "vs1_list")

# loop over groups
dmesvs_mvpp_scores = lapply(seq_along(all_groups_info$group_names), function(gn){
  print(paste0("gn == ", gn))
  # get group name:
  file_name  <- paste0(all_groups_info$group_names[gn], "_mout_", all_groups_info$output_dim_standard[gn])
  
  # loop over the energy and variogram scores
  dmesvs_mvpp_scores_tmp = list()
  for(sc in seq_along(dmesvs_mvpp_scores_list)){
    keep_scores = dmesvs_mvpp_scores_list[sc] 

    # load the actual data and get the data.frame, and filter:
    data_env <- load_and_prepare_data(file_name, dmesvs_mvpp_benchmark_name)
    df <- data_env$df
    dfplot <- filter_models(df, dmesvs_mvpp_model_names, keep_scores)$df %>% subset(score %in% keep_scores)
    # rearrange for plotting
    dfplot <- dfplot %>%
      pivot_longer(cols = starts_with("bootstrap_"), names_to = "bootstrap", values_to = "value") %>%
      # recode some things for the plotting
      mutate(group_name = file_name, 
             group_names_pretty = all_groups_info$group_names_pretty[gn],
             score = ifelse(dmmvsimssh_scores_list[sc] == "es_list", "Energy Score", "Variogram Score"))
    # get prettier names for the models
    new_model_tmp = unlist(lapply(as.character(dfplot$model), change_names))
    dfplot$model = new_model_tmp
    # get prettier names for the bm
    new_bm_tmp = unlist(lapply(as.character(dfplot$benchmark), change_names))
    dfplot$benchmark = new_bm_tmp
    
    # save for this score
    dmesvs_mvpp_scores_tmp[[sc]] = dfplot 
    
  }
  return(bind_rows(dmesvs_mvpp_scores_tmp))
  
  
}) %>% bind_rows()

# get prettier group names for the plot
dmesvs_mvpp_scores$group_names_pretty = factor(gsub("_mout_[0-9][0-9]", "", dmesvs_mvpp_scores$group_names_pretty), levels = all_groups_info$group_names_pretty)
# get prettier model names for the plot
dmesvs_mvpp_scores$model = factor(dmesvs_mvpp_scores$model, levels = c("SSh", "SimSSh", "ECC", "GCA"))

pdf(file = paste0(config$cobase_dir, config$results_folder, "Figures/flos_dmesvs_mvpp.pdf"), width = 12, height = 8)
print(
ggplot(dmesvs_mvpp_scores, aes(model, value)) +
  facet_wrap(~score+group_names_pretty, nrow = 2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = qnorm(alpha), ymax = qnorm(1 - alpha)),
            fill = "gray75", color = "gray75", alpha = alpha) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray25") +
  theme_bw() +
  xlab("Reshuffling method") + ylab("DM test statistic") +
  ggtitle(label = paste0("MVPP method v ", "COBASE-GCA")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(size=16),
        plot.margin = margin(t = 10, r = 33, b = 10, l = 45))
)
dev.off()



########################################
########## dmesvs_parcop (Figure 4):
########################################
# 
print(paste0("working on: dmesvs_parcop"))

dmesvs_parcop_model_names     = c(config$dmesvs_parcop_models)
dmesvs_parcop_benchmark_name  = config$dmesvs_parcop_bm 
dmesvs_parcop_scores_list = c("es_list", "vs1_list")


# loop over groups
dmesvs_parcop_scores = lapply(seq_along(all_groups_info$group_names), function(gn){
  print(paste0("gn == ", gn))
  # get group name:
  file_name  <- paste0(all_groups_info$group_names[gn], "_mout_", all_groups_info$output_dim_standard[gn])

  # loop over the energy and variogram scores
  dmesvs_parcop_scores_tmp = list()
  for(sc in seq_along(dmesvs_parcop_scores_list)){
    keep_scores = dmesvs_parcop_scores_list[sc] 

    # each model has it's own bm, so just loop over the bms
    df2_bb = list()
    for(bb in seq_along(dmesvs_parcop_benchmark_name)){
      # load the actual data and get the data.frame, and filter to keep models that use this bm:
      data_env <- load_and_prepare_data(file_name, dmesvs_parcop_benchmark_name[bb])
      df = data_env$df
      # now filter to keep only the model we need and the score we want:
      df2_bb[[bb]] <- filter_models(df, dmesvs_parcop_model_names[bb], keep_scores)$df
    }
    df2 = bind_rows(df2_bb)
    # rearrange for plotting
    dfplot <- subset(df2, score %in% keep_scores)
    dfplot <- dfplot %>%
      pivot_longer(cols = starts_with("bootstrap_"), names_to = "bootstrap", values_to = "value") %>%
      # recode some things for the plotting
      mutate(group_name = file_name, 
             group_names_pretty = all_groups_info$group_names_pretty[gn],
             score = ifelse(dmesvs_parcop_scores_list[sc] == "es_list", "Energy Score", "Variogram Score"))
    
    # get prettier names for the models
    new_model_tmp = unlist(lapply(as.character(dfplot$model), change_names))
    dfplot$model = new_model_tmp
    # get prettier names for the bm
    new_bm_tmp = unlist(lapply(as.character(dfplot$benchmark), change_names))
    dfplot$benchmark = new_bm_tmp
    
    # save for this score
    dmesvs_parcop_scores_tmp[[sc]] = dfplot 
    
  }
  return(bind_rows(dmesvs_parcop_scores_tmp) )
  
  
}) %>% bind_rows()

# get prettier group names for the plot
dmesvs_parcop_scores$group_names_pretty = factor(gsub("_mout_[0-9][0-9]", "", dmesvs_parcop_scores$group_names_pretty), levels = all_groups_info$group_names_pretty)
# get prettier model names for the plot
dmesvs_parcop_scores$model = factor(dmesvs_parcop_scores$model, levels = unlist(lapply(dmesvs_parcop_model_names, change_names)))



pdf(file = paste0(config$cobase_dir, config$results_folder, "Figures/flos_dmesvs_parcop.pdf"), width = 12, height = 8)
print(
ggplot(dmesvs_parcop_scores, aes(model, value)) +
  facet_wrap(~score+group_names_pretty, nrow = 2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = qnorm(alpha), ymax = qnorm(1 - alpha)),
            fill = "gray75", color = "gray75", alpha = alpha) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray25") +
  theme_bw() +
  xlab("Parametric copula method") + ylab("DM test statistic") +
  ggtitle(label = paste0("COBASE-model v ", "model")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(size=16),
        plot.margin = margin(t = 10, r = 33, b = 10, l = 45)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), col = 'grey75', lty = 'dashed')
)
dev.off()





########################################
########################################
########## crps table (Table B.3):
########################################
# 
crps_model_names_list = c("ens", config$crps_scores_models)


# loop over groups
out_crps <- lapply(seq_along(all_groups_info$group_names), function(gn){
  print(paste0("gn == ", gn))
  # get group name
  file_name  <- paste0(all_groups_info$group_names[gn], "_mout_", all_groups_info$output_dim_standard[gn])
  
  
  # load scores:
  res <- new.env()
  load_score_env(file_name, res)
  
  # keep score data that we want:
  # loop over the models that we want to keep:
  out_scores_tmp <- lapply(seq_along(crps_model_names_list), function(mm) {
    # this is a data.frame with stations in columns and days in rows. take the mean over all days per station    
    tmp = colMeans(res[["crps_list"]][[crps_model_names_list[mm]]])
    # make a data.frame for the table
    tmp2 = data.frame(crps = c(tmp),
                      st_num = names(tmp),
                      model = crps_model_names_list[mm],
                      score = "crps_list") %>% 
      separate(st_num, into = c("var", "st_num"), sep = "-") 
       return(tmp2)
      }) %>% bind_rows()
    
    # get prettier models for plotting
    new_model_tmp = unlist(lapply(as.character(out_scores_tmp$model), change_names))
    out_scores_tmp$model = new_model_tmp
    return(out_scores_tmp )
    
}) %>% bind_rows()

# this is just for checking DPT
# out_crps_df_dpt = out_crps %>% 
#   filter(var %in% c("T_DEWP_10")) %>%
#   merge(., station_info) %>%
#   mutate(Station = st_name) %>%
#   dplyr::select(Station, model, crps)
# 
# ggplot(out_crps_df_dpt, aes(x = model, y = crps, color = Station)) + geom_point() + ggtitle("crps scores")


# this is for making the T2m table
out_crps_df = out_crps %>% 
  # select "obs" which is T2m from AL, or T_DRYB_10 which is T2m from EC
  filter(var %in% c("obs", "T_DRYB_10")) %>%
  # put with the prettier station information
  merge(., station_info) %>%
  # rename and select things for the table
  mutate(Station = st_name) %>%
  dplyr::select(Station, model, crps)

# plot to check 
# ggplot(out_crps_df, aes(x = model, y = crps, color = Station)) + geom_point()


# reformat for the table
out_crps_df = out_crps_df  %>% 
  # Sonnblick is in there twice. remove one.
  distinct(Station, model, .keep_all = TRUE) %>%
  pivot_wider(names_from = model,values_from = crps)

# reorder models for the table
out_crps_df <- out_crps_df[, c("Station", unlist(lapply(as.character(crps_model_names_list[-1]), change_names)))]

# make an xtable object
out_crps_df_tab <- xtable(out_crps_df, digits = 4)
print(out_crps_df_tab, include.rownames = FALSE, 
      file = paste0(config$cobase_dir, config$results_folder, "Figures/crps_table.tex"))

  
########################################
########## esvs table (Table B4):
########################################
# 
# get a list of model names
esvs_model_names_list = c("ens", config$mvpp_scores_models)
# pretty model names in the specific order that we want:
esvs_model_names_list_order = c("SimSSh-R", "SimSSh", "SSh", "ECC", "GCA", "Clayton", "Frank", "Gumbel", "COBASE-GCA", "COBASE-Clayton", "COBASE-Frank", "COBASE-Gumbel" )
mvscores_list = c("es_list", "vs1_list")


# loop over groups:
out_esvs <- lapply(seq_along(all_groups_info$group_names), function(gn){
  print(paste0("gn == ", gn))
  # get group name:
  file_name  <- paste0(all_groups_info$group_names[gn], "_mout_", all_groups_info$output_dim_standard[gn])
  
  
  # load scores:
  res <- new.env()
  load_score_env(file_name, res)
  
  # loop over the models that we want to keep:
  out_scores_tmp <- lapply(seq_along(esvs_model_names_list), function(mm) {
    # loop over the scores:
    lapply(mvscores_list, function(s){
      # each item is a list of the score for each group. take a mean of it
      tmp = mean(res[[s]][[esvs_model_names_list[mm]]])
      # make a data.frame for the table
      tmp2 = data.frame(value = c(tmp),
                        group = all_groups_info$group_names[gn],
                        group_pretty = all_groups_info$group_names_pretty[gn],
                        model = esvs_model_names_list[mm],
                        score = ifelse(s == "es_list", "Energy Score", "Variogram Score")) 
      return(tmp2)
    }) %>% bind_rows()
    
  }) %>% bind_rows() 
  
  # get prettier models for plotting
  new_model_tmp = unlist(lapply(as.character(out_scores_tmp$model), change_names))
  out_scores_tmp$model = new_model_tmp
  return(out_scores_tmp )
  
}) %>% bind_rows() %>% unique()


# make a plot to check:
# ggplot(out_esvs, aes(x = model, y = value, color = model)) +
#   geom_point() +
#   facet_wrap(~score+group, nrow = 2) +
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
#         axis.text.y = element_text(size = 12),
#         text=element_text(size=16),
#         plot.margin = margin(t = 10, r = 33, b = 10, l = 45))



# reformt for the table:
out_esvs_formatted <- out_esvs %>%
  # rename scores:
  mutate(score_abbrev = ifelse(score == "Energy Score", "ESS", "VSS")) %>%
  # get prettier group names
  mutate(group_score = paste(group_pretty, score_abbrev, sep = "_")) %>%
  # remove variables we don't need in the table
  dplyr::select(-score, -score_abbrev, -group, -group_pretty) %>%
  pivot_wider(names_from = group_score, values_from = value)

# make into an xtable object
out_esvs_formatted_tab <- xtable(out_esvs_formatted, digits=4)
## save the scores:
print(out_esvs_formatted_tab, include.rownames = FALSE, 
      file = paste0(config$cobase_dir, config$results_folder, "Figures/esvs_table.tex"))


