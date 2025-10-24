source(paste0(config$cobase_dir, "/Postprocessing/Utilities/load_data_util.R"))
library(tidyr)
library(xtable)

save_table_full <- function(file_name, model_names, scores, save_folder, show_sd = TRUE)
{

    # Retrieve the scores
    res <- new.env()
    load_score_env(file_name, res)

    table_preamble <- sprintf("
    \\begin{table}[!ht]
    \\centering
    \\makebox[\\textwidth]{
    \\begin{tabular}{|ll%s|}
    \\hline
    \\multicolumn{2}{|l|}{}  %s\\\\\\hline\n", 
    strrep("|l", length(model_names)),
    paste(paste0(sapply(model_names, function(x) paste0(" & \\textbf{", x, "}"))), collapse = ""))

    score_names <- vapply(scores, pretty_score_name, character(1))

    for (i in 1:length(scores))
    {
    score <- scores[i]
    score_name <- score_names[i]

    scores_per_model <- res[[score]]
    
    ref_model <- model_names[2]
    
    if (length(dim(scores_per_model[[ref_model]])) == 2)
    { # Average over the first dimension
        d <- dim(scores_per_model[[ref_model]])[2]
        
        table_preamble <- paste0(table_preamble, "\\multicolumn{1}{|c|}{\\multirow{", toString(d), "}{*}{", score_name, "}}\n")
        
        for (dd in 1:d)
        {
        if (dd != 1)
        {
            table_preamble <- paste0(table_preamble, "\\multicolumn{1}{|l|}{} \n")
        }
        table_preamble <- paste0(table_preamble, " & ", toString(dd))
        for (model in model_names)
        {
            val <- mean(scores_per_model[[model]][,dd])
            sd_val <- sd(scores_per_model[[model]][,dd])
            
            if (show_sd) 
                sd_str <- "\\pm " + toString(round(sd_val, 2)) + "$"
            else
                sd_str <- "$"
            
            table_preamble <- paste0(table_preamble, "& $", toString(round(val, 2)), sd_str)
        }
        if (dd != d)
        {
            table_preamble <- paste0(table_preamble, "\\\\ \\cline{2-", toString(length(model_names) + 2), "}\n")
        } else
        {
            table_preamble <- paste0(table_preamble, "\\\\ \\hline\n")
        }
        
        }
        
    } else 
    {
        table_preamble <- paste0(table_preamble, "\\multicolumn{2}{|c|}{", score_name, "}\n")
        
        for (model in model_names)
        {
        val <- mean(scores_per_model[[model]])
        sd_val <- sd(scores_per_model[[model]])

        if (show_sd) 
            sd_str <- "\\pm " + toString(round(sd_val, 2)) + "$"
        else
            sd_str <- "$"

        table_preamble <- paste0(table_preamble, " & $", toString(round(val, 2)), sd_str)
        }
        table_preamble <- paste0(table_preamble, "\\\\\\hline \n")
        
    }
    }

    table_preamble <- paste0(table_preamble, sprintf("\\end{tabular}
    }
    \\end{table}"))

    saveFile <- paste0(save_folder, "table_", file_name, ".txt")

    sink(saveFile)
    cat(table_preamble)
    sink()

}

