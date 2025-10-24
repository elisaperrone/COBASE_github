# creating all results with the config file:

# install packages: config, scoringRules, remotes, verification, car, ggpubr 


# 1. Load the config file:
#setwd("./COBASE_github")
config <- config::get()

# 2. Univariate PostProcessing
# calls uvpp_util.R to do the univariate post-processing and calcualte the CRPS
# calss sim_util.R to calculate the sim_matrix
# saves output
source(paste0(config$cobase_dir, "/Postprocessing/univariate.R"))

# 3. post-process the NWP forecasts. Verify the post-processed forecasts. Calculate DM statistics
# calling postprocess.R does the following:
# - runs the function "postprocess_all" from mvpp.R
# - in "postprocess_all" we define a function "mvpp_adjusted" that actually post-processes the methods: 
#   - emos.q, emos.r
#   - ecc.q
#   - ssh.i, ssh.sim (emos.q), ssh.sim (emos.r)
#   - gca, COBASE-gca (emos.q)
#   - clayton, gumbel, frank, COBASE-clayton (emos.q), COBASE-gumbel (emos.q), COBASE-frank (emos.q)
# runs the function "compute_DM_scores" from DM_util.R - computes scores and saves them

source(paste0(config$cobase_dir, "/Postprocessing/multivariate.R"))


# 3. creating plots:
# calling create_plots_paper.R does the following:
# - defines a function "create_all_results"
# - "create_all_results" calls the function "create_dm_plots_single_benchmark" from dm_plots.R
source(paste0(config$cobase_dir, "/Creating Plots/create_plots_paper_selection.R"))






