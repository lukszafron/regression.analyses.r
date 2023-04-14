#! /usr/bin/env Rscript

cat("Program path:", unlist(strsplit(grep(commandArgs(), pattern = "file=", value = T), split = "="))[2], "\n")

args <- commandArgs(trailingOnly = T)

if(length(args) != 14) {cat(paste("The number of provided arguments is incorrect:", 
                                      length(args), "instead of 14.
                                       The arguments should be placed in the following order:
                                       1 - a path to a semicolon-separated CSV file containing the variables to be analyzed,
                                       2 - a variable in the aforementioned CSV file where sample names are stored,
                                       3 - are the veriables to be analyzed continuous ('TRUE', 'FALSE'),
                                       4 - a path to a semicolon-separated CSV file containing clinico-pathological and follow-up data,
                                       5 - a variable in the aforementioned CSV file where sample names are stored,
                                       6 - a vector listing any number of dependent clinico-pathological/follow-up variables to be used in the following format: 'discrete_var1;discrete_var2;continuous_time_var1:discrete_censoring_var1;continuous_time_var2:discrete_censoring_var2',
                                       7 - a vector listing any number of independent discrete variables to be used in the following format: 'independent_discrete_var1;independent_discrete_var2', 'NA' if none,
                                       8 - a vector listing any number of independent continuous variables to be used in the following format: 'independent_continuous_var1;independent_continuous_var2', 'NA' if none,
                                       9 - a vector listing one or two grouping variables to be used in the following format: 'grouping_var1;grouping_var2', 'NA' if none,
                                       10 - a number of threads to be used for parallel computations,
                                       11 - a vector listing the types of regression analyses to be performed, 'coxph' for Cox regression, 'lrm' for logistic regression, separated with a semicolon,
                                       12 - a path to the working directory,
                                       13 - a statistical significance level (alpha) to be used,
                                       14 - should the Benjamini-Hochberg correction for multiple comparisons be used ('TRUE', 'FALSE').\n"))
  stop("Incorrect number of arguments.")}

args

library(survival)
library(rms)
library(foreach)
library(sjmisc)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(riskRegression)
library(ggplot2)
library(ggtext)
library(pdftools)
library(pals)
library(grid)
library(gridExtra)
library(survminer)
library(stringi)
library(tidyr)
library(epiDisplay)
library(data.table)
library(doMC)
threads <- as.numeric(args[10])
registerDoMC(threads)

main.table <- as.data.frame(fread(file = args[1], header = T, check.names = T, sep = ";", na.strings = c("", NA)))
main.table <- main.table %>% transmute(across(.cols = everything(), .fns = function(x){if(is.character(x)) {gsub(x, pattern = "(\\.| )+", replacement = "_")} else {x}}))

main.table.path <- unlist(strsplit(args[1], split = "/"))
main.table.name <- main.table.path[length(main.table.path)]
main.samples <- args[2]
main.wo.samples <- main.table %>% dplyr::select(-main.samples)
main.list <- colnames(main.table)[colnames(main.table) != main.samples]
is.continuous <- as.logical(args[3])
if(!is.continuous) {
  for(i in colnames(main.table)) {
  main.table[,i] <- as.factor(main.table[,i])}} else if(!all(sapply(main.wo.samples, is.numeric))) {
  non.numeric <- names(main.wo.samples)[!sapply(main.wo.samples, is.numeric)]
  warning(paste0("The main table contains some non-numeric values which are not allowed for continuous variables:\n", paste(non.numeric, collapse = ", "),".\n These variables will be deleted."))
  main.table <- main.table %>% dplyr::select(-non.numeric)
  }

df <- as.data.frame(fread(file = args[4], header = T, check.names = T, sep = ";", na.strings = c("", NA)))
df <- df %>% transmute(across(.cols = everything(), .fns = function(x){if(is.character(x)) {gsub(x, pattern = "(\\.| )+", replacement = "_")} else {x}}))
df.table.path <- unlist(strsplit(args[4], split = "/"))
df.table.name <- df.table.path[length(df.table.path)]
df.samples <- args[5]
df <- df %>% dplyr::select(-colnames(main.table)[colnames(main.table) %in% colnames(df) & colnames(main.table) %nin% main.samples])
df <- merge(x = main.table, y = df, by.x = main.samples, by.y = df.samples)
df <- df %>% transmute(across(.fns = function(x) {if(class(x) == "character") {as.factor(x)} else {x}}))

dep.vars.vec <- args[6]
other.discrete.factors.vec <- args[7]
other.discrete.factors <- unlist(strsplit(other.discrete.factors.vec, split = ";"))
if(all(other.discrete.factors != "NA")) {for(i in other.discrete.factors) {df[[i]] <- as.factor(df[[i]])}}
other.continuous.factors.vec <- args[8]
other.continuous.factors <- unlist(strsplit(other.continuous.factors.vec, split = ";"))
other.factors <- c(other.discrete.factors, other.continuous.factors)
other.factors <- other.factors[other.factors != "NA"]
if(length(other.factors) == 0) {variate <- "univar"} else {variate <- "multivar"}

dep.vars <- unlist(strsplit(dep.vars.vec, split = ";"))
sub.vars.vec <- args[9]
sub.vars <- unlist(strsplit(sub.vars.vec, split = ";"))

atypes.vec <- args[11]
atypes <- unlist(strsplit(atypes.vec, split = ";"))
workdir <- args[12]
p.value <- as.numeric(args[13])
fdr <- as.logical(args[14])

dir.create(path = workdir, recursive = T, showWarnings = F)
setwd(workdir)
unlink(list.files(pattern = "^(riskRegression\\.)([0-9]+)(.*\\.pdf\\.tmp$)"))

if(length(sub.vars) >= 1 & all(sub.vars != "ALL_SAMPLES")) {
for(i in sub.vars) {df[[i]] <- as.factor(df[[i]])}
sub.vars.list <- foreach(sub.var = sub.vars) %do% {
  grep(levels(df[[sub.var]]), pattern = "^$", invert = T, value = T)
}
names(sub.vars.list) <- sub.vars} else if(length(sub.vars) == 1 & sub.vars == "ALL_SAMPLES") {
sub.vars.list <- list("ALL_SAMPLES" = "ALL_SAMPLES")} else 
  {stop("The subsetting vars list seems to be incorrect.")}

file.name.prefix <- paste("Reg_anal", paste(atypes, collapse = ","), variate, main.table.name, paste("cont", is.continuous, sep = ":"), df.table.name, paste0("dep_vars:", paste(dep.vars, collapse = ",")), paste0("group_vars:", paste(sub.vars, collapse = ",")), paste0("other_vars:", paste(other.factors, collapse = "+")), sep = ".")
file.name.prefix <- gsub(file.name.prefix, pattern = "BAM_files_with_duplicates", replacement = "BAM_w.dups")
file.name.prefix <- gsub(file.name.prefix, pattern = "BAM_files_without_duplicates", replacement = "BAM_wo.dups")
file.name.prefix <- gsub(file.name.prefix, pattern = "TRUE", replacement = "T")
file.name.prefix <- gsub(file.name.prefix, pattern = "FALSE", replacement = "F")

if(length(other.factors) == 0) {file.name.prefix <- paste0(file.name.prefix, "NA")}
file.name.prefix <- sub(stri_reverse(substr(stri_reverse(file.name.prefix), start = 1, stop = 212)), pattern = "^\\.*", replacement = "")

# save.image(paste(file.name.prefix, "RData", sep = "."))

if(! file.exists(paste(file.name.prefix, "res", "RData", sep = "."))) {
  unlink("Regression_analyses_progress.log")
    
  res <- as.data.frame(rbindlist(foreach(index = seq(main.list)) %dopar% {
    gene <- main.list[index]
    if(index == 1 | (index %% 100) == 0 | index == length(main.list)) {
      cat("Performing the regression analysis for the entry no. ", index, " out of ",
          length(main.list), " entries ", "(",
          sprintf(fmt = "%0.1f%%", index/length(main.list)*100), "), date: ", as.character(Sys.time()), "...\n",
          sep = "", file = "Regression_analyses_progress.log", append = T)}
    
    test.summary <- function(test) {
      if(all(class(test) != "character" & class(test) != "NULL")) {
        if(length(test$coefficients) > 0) {
          test$call$formula <- deparse(eval(parse(text = test$call$formula)))
          test$call$formula <- gsub(paste(test$call$formula, collapse = ""), pattern = " {2,}", replacement = " ")
          if(exists("sub.value")) {test$call$data <- sub.value} else {test$call$data <- "full_table"}
          
          if(is.null(test$conf.int)) {
            if(!is.null(lrm1)) {
              if(!lrm1$fail) {
                res.odds <- as.matrix(summary(lrm1))
                res.odds <- exp(res.odds[seq(1,nrow(res.odds),2),, drop = F])
                gene.names <- sub(rownames(res.odds), pattern = " - ", replacement = "")
                gene.names <- foreach(gname = gene.names, .combine = c) %do% {
                  if(gname != gene) {sub(gname, pattern = ":[^:]+$", replacement = "")} else {gname}}
                
                rownames(res.odds) <- gene.names
                test$conf.int <- res.odds
                if(all(str_contains(rownames(test$conf.int), pattern = rownames(test$coefficients)))) {
                  row.order <- foreach(row = rownames(test$coefficients), .combine = c) %do% {
                    grep.val <- grep(rownames(test$conf.int), pattern = paste0("^", row, "$"))
                    if(length(grep.val) != 1) {stop("Row names are missing or duplicated.")}
                    grep.val}
                  if(length(row.order) != length(rownames(test$conf.int))) {stop("Rows reordering failed.")}
                  test$conf.int <- test$conf.int[row.order,, drop = F]
                } else {stop("Independent variable names do not match.")}
              } else {
                glm.odds <- tryCatch(expr = logistic.display(glm0), error = function(e) {return(NULL)})
                if(!is.null(glm.odds)) {
                  if(nrow(glm.odds$table) > 2) {
                    glm.odds <- logistic.display(glm0, simplified = T)
                    test$conf.int <- glm.odds$table
                  } else {
                    glm.odds$table <- as.data.frame(glm.odds$table)
                    glm.odds.rowname <- sub(rownames(glm.odds$table)[1], pattern = " \\(cont. var.\\)", replacement = "")
                    if(length(unlist(strsplit(glm.odds.rowname, split = " "))) >1) {
                      glm.odds.rowname <- paste(sub(unlist(strsplit(glm.odds.rowname, split = " "))[1:2], pattern = ":$", replacement = ""), collapse = "")
                    }
                    glm.odds.rowvalues <- t(as.matrix(as.numeric(unlist(strsplit(gsub(gsub(as.character(glm.odds$table$`OR(95%CI)`), pattern = "[()]", replacement = ","), pattern = "[ ]", replacement = "")[1], split = ",")))))
                    rownames(glm.odds.rowvalues) <- glm.odds.rowname
                    colnames(glm.odds.rowvalues) <- c("OR", "lower95ci", "upper95ci")
                    test$conf.int <- glm.odds.rowvalues
                  }} else {
                    rows.n <- length(glm0$coefficients[names(glm0$coefficients) != "(Intercept)"])
                    if(rows.n > 1) {
                      glm.odds <- foreach(n = seq(1,rows.n), .combine = rbind) %do% {c(NA, NA, NA)}} else {
                        glm.odds <- t(as.data.frame(rep(NA, 3)))}
                    rownames(glm.odds) <- names(glm0$coefficients)[names(glm0$coefficients) != "(Intercept)"]
                    test$conf.int <- glm.odds
                  }
              }} else {
                glm.odds <- tryCatch(expr = logistic.display(glm0), error = function(e) {return(NULL)})
                if(!is.null(glm.odds)) {
                  if(nrow(glm.odds$table) > 2) {
                    glm.odds <- logistic.display(glm0, simplified = T)
                    test$conf.int <- glm.odds$table
                  } else {
                    glm.odds$table <- as.data.frame(glm.odds$table)
                    glm.odds.rowname <- sub(rownames(glm.odds$table)[1], pattern = " \\(cont. var.\\)", replacement = "")
                    if(length(unlist(strsplit(glm.odds.rowname, split = " "))) >1) {
                      glm.odds.rowname <- paste(sub(unlist(strsplit(glm.odds.rowname, split = " "))[1:2], pattern = ":$", replacement = ""), collapse = "")
                    }
                    glm.odds.rowvalues <- t(as.matrix(as.numeric(unlist(strsplit(gsub(gsub(as.character(glm.odds$table$`OR(95%CI)`), pattern = "[()]", replacement = ","), pattern = "[ ]", replacement = "")[1], split = ",")))))
                    rownames(glm.odds.rowvalues) <- glm.odds.rowname
                    colnames(glm.odds.rowvalues) <- c("OR", "lower95ci", "upper95ci")
                    test$conf.int <- glm.odds.rowvalues
                  }} else {
                    rows.n <- length(glm0$coefficients[names(glm0$coefficients) != "(Intercept)"])
                    if(rows.n > 1) {
                      glm.odds <- foreach(n = seq(1,rows.n), .combine = rbind) %do% {c(NA, NA, NA)}} else {
                        glm.odds <- t(as.data.frame(rep(NA, 3)))}
                    rownames(glm.odds) <- names(glm0$coefficients)[names(glm0$coefficients) != "(Intercept)"]
                    test$conf.int <- glm.odds
                  }}
          }
          rownames(test$coefficients) <- rows.rename(test$coefficients)
          rownames(test$conf.int) <- rows.rename(test$conf.int)
          df1 <- t(as.data.frame(c(test$call$formula, test$call$data)))
          rownames(df1) <- NULL
          colnames(df1) <- c("Formula", "Data")
          if(atype == "coxph") {
            test.res <- merge(x = test$conf.int[,c(1,3:4), drop = F], y = setNames(as.data.frame(test$coefficients[,5, drop = F]), nm = "p-values"), by.x = 0, by.y = 0, sort = F)} else
              if(atype == "lrm") {
                if(class(lrm1)[1] == "lrm") {
                  if(!lrm1$fail) {
                    test.res <- merge(x = test$conf.int[,c(4,6:7), drop = F], y = setNames(as.data.frame(test$coefficients[,4, drop = F]), nm = "p-values"), by.x = 0, by.y = 0, sort = F)} else
                      if(lrm1$fail & nrow(test$conf.int) >= 1) {test.res <- merge(x = test$conf.int[,c(1:3), drop = F], y = setNames(as.data.frame(test$coefficients[,4, drop = F]), nm = "p-values"), by.x = 0, by.y = 0, sort = F)}} else
                        if(is.null(lrm1) & nrow(test$conf.int) >= 1){test.res <- merge(x = test$conf.int[,c(1:3), drop = F], y = setNames(as.data.frame(test$coefficients[,4, drop = F]), nm = "p-values"), by.x = 0, by.y = 0, sort = F)} else
                        {stop("The conf. int. table contains no rows.")}}
          
          colnames(test.res)[1] <- "Factor"
          if(any(is.na(test.res[,3:4]))) {
            test.res <- merge(x = test.res, y = gsub(as.matrix(test.res[,3:4] %>% unite(col = "95% CI", sep = "-")), pattern = "(.*)", replacement = "[\\1]"), by.x = 0, by.y = 0)} else {
              test.res <- merge(x = test.res, y = gsub(as.matrix(round(test.res[,3:4],3) %>% unite(col = "95% CI", sep = "-")), pattern = "(.*)", replacement = "[\\1]"), by.x = 0, by.y = 0)}
          test.res <- test.res[,c(2,3,7,6,4,5)]
          colnames(test.res) <- c("Factor", "HR/OR", "95% CI", "p-value", "Lower 0.95 CI", "Upper 0.95 CI")
          if(atype == "coxph") {
            test.res[["N"]] <- test[["n"]]
            test.res[["Events no."]] <- test[["nevent"]] } else
              if(atype == "lrm") {
                test.res[["N"]] <- test$df.null + 1
                test.res[["Events no."]] <- sum(glm0$model[[dep.var]] == levels(glm0$model[[dep.var]])[2])
              }} else {
                df1 <- as.data.frame(t(c(gsub(paste(deparse(final.form), collapse = ""), pattern = " {2,}", replacement = " "), sub.value)))
                rownames(df1) <- NULL
                colnames(df1) <- c("Formula", "Data")
                if(length(test$coefficients) == 0) {test.res <- as.data.frame(t(c(paste0("#", gene), rep("All independent vars are constant.",2), 69, rep("All independent vars are constant.",4))))} else 
                {
                  stop("The regression analysis failed due to an unknown error.")
                }
                colnames(test.res) <- c("Factor", "HR/OR", "95% CI", "p-value", "Lower 0.95 CI", "Upper 0.95 CI", "N", "Events no.")
              }} else {
                df1 <- as.data.frame(t(c(gsub(paste(deparse(final.form), collapse = ""), pattern = " {2,}", replacement = " "), sub.value)))
                rownames(df1) <- NULL
                colnames(df1) <- c("Formula", "Data")
                
                if(is.null(test)) {
                  test.res <- as.data.frame(t(c(paste0("#", gene), rep("Too few obs.",2), 69, rep("Too few obs.",4))))} else
                    if(test[1] == "ERROR") {
                      test.res <- as.data.frame(t(c(paste0("#", gene), rep("Too many vars, constant vars or vars collinearity",2), 69, rep("Too many vars, constant vars or vars collinearity",4))))} else
                        if(length(test$coefficients) == 0) {
                          test.res <- as.data.frame(t(c(paste0("#", gene), rep("All independent vars are constant.",2), 69, rep("All independent vars are constant.",4))))} else 
                        {
                          stop("The regression analysis failed due to an unknown error.")
                        }
                colnames(test.res) <- c("Factor", "HR/OR", "95% CI", "p-value", "Lower 0.95 CI", "Upper 0.95 CI", "N", "Events no.")
              }
      df.full <- cbind(df1, test.res)
      suppressWarnings(rm(lrm1, dd, glm0))
      #    options(datadist = NULL)
      return(df.full)}
    rows.rename <- function(x) {
      new.rownames <- foreach(i = rownames(x), .combine = c) %do% {
        if(str_contains(x = i, pattern = other.factors, logic = "OR") & i %nin% other.factors & i %nin% gene) {
          suffix <- unlist(base::strsplit(i, split = other.factors[str_contains(x = i, pattern = other.factors)]))[2]
          varname <- sub(i, pattern = paste0(suffix, "$"), replacement = "")
          full.name <- paste(varname, suffix, sep = ":")
          full.name } else if(i %in% other.factors) {i} else if(i == gene) {
            if(is.continuous){
              paste0("#", gene)} else {stop(paste("The", i, "gene name does not contain the category suffix."))}
          } else if(str_contains(i, pattern = gene)) {
            suffix <- unlist(strsplit(i, split = gene))[2]
            varname <- sub(i, pattern = paste0(suffix, "$"), replacement = "")
            full.name <- paste(paste0("#", varname), suffix, sep = ":")
            full.name
          } else {stop(paste("The", i, "factor seems invalid."))}}
      return(new.rownames)
    }
    lrm.prep <- function(data) {
      lrm1 <<- tryCatch(expr = lrm(final.form, data = data), error = function(e) {return(NULL)})
      glm0 <<- tryCatch(expr = glm(final.form, data = data, family = "binomial"), error = function(e) {return(NULL)})
      if(!is.null(glm0)) {
        glm1 <- summary(glm0)
        glm1$coefficients <- glm1$coefficients[rownames(glm1$coefficients) != "(Intercept)", , drop = F]
        test.summary(test = glm1)
      } else {
        test.summary(test = "ERROR")}}
    
    foreach(atype = atypes, .combine = rbind) %do% {
      foreach(dep.var = dep.vars, .combine = rbind) %do% {
        dependent.var <- unlist(strsplit(dep.var, split = ":"))[1]
        censoring.var <- unlist(strsplit(dep.var, split = ":"))[2]
        
        if(!is.na(censoring.var)) {
          df[[censoring.var]] <- as.factor(df[[censoring.var]])
          df.sub <- df[!is.na(df[[censoring.var]]), , drop = F]} else
            if(is.na(censoring.var)) {
              df[[dependent.var]] <- as.factor(df[[dependent.var]])
              df.sub <- df
            }
        df.sub <- df.sub[!is.na(df.sub[[dependent.var]]), , drop = F]
        selected.columns <- c(gene, dependent.var, censoring.var, other.factors, sub.vars)
        selected.columns <- selected.columns[!is.na(selected.columns)]
        selected.columns <- selected.columns[selected.columns != "ALL_SAMPLES"]
        df.sub <- df.sub %>% dplyr::select(all_of(selected.columns))

        if(length(levels(df.sub[[censoring.var]])) %nin% c(1,2) & !is.na(censoring.var)) {
          stop(paste("The censoring variable:", censoring.var, "is invalid since it has 0 or more than 2 levels."))}
        
        if(length(levels(df.sub[[dependent.var]])) %nin% c(1,2) & is.na(censoring.var)) {
          stop(paste("The dependent variable:", dependent.var, "is invalid since it has 0 or more than 2 levels."))}
        
        if(atype == "coxph" & !is.na(censoring.var)) {
          Surv.form <- paste0("Surv(", dependent.var, ",", censoring.var, " == '", levels(df.sub[[censoring.var]])[2], "')")
          Factors.form <- sub(paste(gene, paste(other.factors, collapse = " + "), sep = " + "), pattern = " \\+ $", replacement = "")
          final.form <- as.formula(paste(Surv.form, Factors.form, sep = " ~ "))
          
          if(length(sub.vars.list) == 1 & all(names(sub.vars.list) == "ALL_SAMPLES")) {
            coxm <- tryCatch(expr = summary(with(df.sub,coxph(formula = final.form, df.sub))), error = function(e){return(NULL)})
            sub.value <- "full_table"
            test.summary(test = coxm)} else if(length(sub.vars.list) == 1) {
              subs <- foreach(sub1 = sub.vars.list[[1]], .combine = c) %do% {paste(names(sub.vars.list[1]), paste0("'", sub1, "'"), sep = "==")}
              subs <- c("full_table", subs)
              foreach(sub.value = subs, .combine = rbind) %do% {if(sub.value == "full_table") {
                coxm <- tryCatch(expr = summary(with(df.sub,coxph(formula = final.form, df.sub))), error = function(e){return(NULL)})
              } else {coxm <- tryCatch(expr = summary(with(df.sub,coxph(formula = final.form, subset(df.sub, subset = eval(parse(text = sub.value)))))), error = function(e){return(NULL)})}
                test.summary(test = coxm)}} else if(length(sub.vars.list) == 2) {
                  subs <- unique(foreach(sub1 = sub.vars.list[[1]], .combine = c) %do% {foreach(sub2 = sub.vars.list[[2]], .combine = c) %do% {c(paste(names(sub.vars.list[1]), paste0("'", sub1, "'"), sep = "=="), paste(names(sub.vars.list[2]), paste0("'", sub2, "'"), sep = "=="), paste(paste(names(sub.vars.list[1]), paste0("'", sub1, "'"), sep = "=="), paste(names(sub.vars.list[2]), paste0("'", sub2, "'"), sep = "=="), sep = " & "))}})
                  subs <- c("full_table", subs)
                  foreach(sub.value = subs, .combine = rbind) %do% {if(sub.value == "full_table") {
                    coxm <- tryCatch(expr = summary(with(df.sub,coxph(formula = final.form, df.sub))), error = function(e){return(NULL)})
                  } else {coxm <- tryCatch(expr = summary(with(df.sub,coxph(formula = final.form, subset(df.sub, subset = eval(parse(text = sub.value)))))), error = function(e){return(NULL)})}
                    test.summary(test = coxm)}
                } else {stop("The number of grouping variables is incorrect.")}
          
        } else if(atype == "lrm" & is.na(censoring.var)){
          dd <<- datadist(df.sub, adjto.cat = "first")
          options(datadist="dd")
          Factors.form <- sub(paste(gene, paste(other.factors, collapse = " + "), sep = " + "), pattern = " \\+ $", replacement = "")
          final.form <- as.formula(paste(dependent.var, Factors.form, sep = " ~ "))
          
          if(length(sub.vars.list) == 1 & all(names(sub.vars.list) == "ALL_SAMPLES")) {
            sub.value <- "full_table"
            lrm.prep(data = df.sub)
          } else if(length(sub.vars.list) == 1) {
            subs <- foreach(sub1 = sub.vars.list[[1]], .combine = c) %do% {paste(names(sub.vars.list[1]), paste0("'", sub1, "'"), sep = "==")}
            subs <- c("full_table", subs)
            foreach(sub.value = subs, .combine = rbind) %do% {if(sub.value == "full_table") {
              lrm.prep(data = df.sub)
            } else {lrm.prep(data = subset(df.sub, subset = eval(parse(text = sub.value))))}}} else
              if(length(sub.vars.list) == 2) {
                subs <- unique(foreach(sub1 = sub.vars.list[[1]], .combine = c) %do% {foreach(sub2 = sub.vars.list[[2]], .combine = c) %do% {c(paste(names(sub.vars.list[1]), paste0("'", sub1, "'"), sep = "=="), paste(names(sub.vars.list[2]), paste0("'", sub2, "'"), sep = "=="), paste(paste(names(sub.vars.list[1]), paste0("'", sub1, "'"), sep = "=="), paste(names(sub.vars.list[2]), paste0("'", sub2, "'"), sep = "=="), sep = " & "))}})
                subs <- c("full_table", subs)
                foreach(sub.value = subs, .combine = rbind) %do% {
                  if(sub.value == "full_table") {
                    lrm.prep(data = df.sub)
                  } else {lrm.prep(data = subset(df.sub, subset = eval(parse(text = sub.value))))}}
              } else {stop("The number of grouping variables is incorrect.")}
        }}}
  }, use.names = T))

res[["p-value"]] <- as.numeric(as.character(res[["p-value"]]))
res[["N"]] <- as.numeric(as.character(res[["N"]]))
res[["Events no."]] <- as.numeric(as.character(res[["Events no."]]))

save(list = "res", file = paste(file.name.prefix, "res", "RData", sep = "."))

} else {
  cat("Loading a pre-existing RData file:", paste(file.name.prefix, "res", "RData...\n", sep = "."))
  load(file = paste(file.name.prefix, "res", "RData", sep = "."))
}

file.name.prefix <- paste(file.name.prefix, paste("pval", p.value, sep = ":"), sep = ".")
file.name.prefix <- paste(file.name.prefix, paste("fdr", fdr, sep = ":"), sep = ".")
file.name.prefix <- gsub(file.name.prefix, pattern = "TRUE", replacement = "T")
file.name.prefix <- gsub(file.name.prefix, pattern = "FALSE", replacement = "F")
file.name.prefix <- sub(stri_reverse(substr(stri_reverse(file.name.prefix), start = 1, stop = 216)), pattern = "^\\.*", replacement = "")

if(fdr) {
res$padj <- p.adjust(res[["p-value"]], method = "fdr")} else {
res$padj <- res[["p-value"]]
}
res <- res[,c(1:6,11,7:10)]
res.sig <- res[grepl(res[["Factor"]], pattern = "^#.*") & res[["padj"]] < p.value & !is.na(res[["padj"]]) & !is.na(res[["HR/OR"]]), , drop = F]
res.united.cols <- res %>% unite(Formula, Data, col = "United.col", sep = "#") %>% dplyr::select("United.col")
res.sig.united.cols <- res.sig %>% unite(Formula, Data, col = "United.col", sep = "#") %>% dplyr::select("United.col")
res.sig.all <- res[res.united.cols[,1] %in% res.sig.united.cols[,1], , drop = F]
res.sig.all <- res.sig.all[res.sig.all[["padj"]] < p.value & !is.na(res.sig.all[["padj"]]) & !is.na(res.sig.all[["HR/OR"]]), , drop = F]
if(is.null(res.sig.all)) {res.sig.all <- setNames(data.frame(t(rep(c(""), 11))), nm = colnames(res))}
res.failed <- subset(res, subset = `p-value` == 69 | is.na(`p-value`))
res.failed <- subset(res.failed, subset = grepl(Factor, pattern = "(^#.*|^Too (many vars|few obs.)|^All independent vars are constant.)"))

dep.var.cats <- sub(dep.vars, pattern = "^.*:", replacement = "")
dep.var.conts <- sub(grep(dep.vars, pattern = ".*:.*", value = T), pattern = ":.*$", replacement = "")
df <- df[df %>% dplyr::select(all_of(c(dep.var.conts))) %>% is.na %>% rowSums() == 0, , drop = F]

if(all(sub.vars != "ALL_SAMPLES")) {
  df <- df %>% mutate(across(.cols = c(dep.var.cats, sub.vars), .fns = as.factor))
  N <- with(df, aggregate(df %>%
                            dplyr::select(dep.var.cats, dep.var.conts, other.factors), 
                          by = sapply(sub.vars, FUN = function(x) {get(x)}, simplify = F), FUN = length))[,dep.var.cats[1]]
  
  summary_f <- function(x) {
    summary.orig <- summary(x)
  if(length(summary.orig) == 7) {
    return(summary.orig)} else
      if(length(summary.orig) == 6) {
        return(c(summary.orig, "NA's" = 0))} else {
          stop("The summary_f function return an unknown error.")}}
  
  df.summary <- with(df, aggregate(df %>%
                                     dplyr::select(dep.var.cats, dep.var.conts, other.factors), 
                                   by = sapply(sub.vars, FUN = function(x) {get(x)}, simplify = F), FUN = function(y) {if(is.factor(y)) {table(y, useNA = "always")} else {summary_f(y)}})) %>%
    as.matrix() %>% as.data.frame(stringsAsFactors = F) %>% mutate(across(.cols = !sub.vars, .fns = as.numeric))

  valid.cols <- foreach(i = colnames(df.summary), .combine = c) %do% {
     if(i %in% sub.vars) {i} else
     if(is.na(sum(df.summary[i]))) {i} else
     if(str_contains(x = i, pattern = sub(dep.var.conts, pattern = "(^.*$)", replacement = "\\1\\."), logic = "OR")) {
       if(!grepl(i, pattern = "\\.NA's$")) {i}} else
     if(sum(df.summary[i]) >0) {i}
      }
  df.summary <- df.summary[valid.cols]
  
  df.summary.fin <- cbind(df.summary, N)
  df.summary.fin <- df.summary.fin[c(1:length(sub.vars), ncol(df.summary.fin), (length(sub.vars)+1):(ncol(df.summary.fin)-1))]
  df.summary.fin <- as.data.frame(t(df.summary.fin[with(df.summary.fin, order(get(sub.vars))),]), stringsAsFactors = F)
  df.summary.fin[c((length(sub.vars)+1):nrow(df.summary.fin)),] <- df.summary.fin[c((length(sub.vars)+1):nrow(df.summary.fin)),] %>% 
    mutate(across(.fns = as.numeric))

  N <- nrow(df)
  df.summary <- df %>% dplyr::select(dep.var.cats, dep.var.conts, other.factors) %>% 
    sapply(., summary) %>% unlist() %>% as.data.frame()
  df.summary.fin.all <- rbind(df.summary, N)
  rownames(df.summary.fin.all)[nrow(df.summary.fin.all)] <- "N"
  df.summary.fin.all <- df.summary.fin.all[c(nrow(df.summary.fin.all), 1:(nrow(df.summary.fin.all)-1)), , drop = F]
  rownames(df.summary.fin) <- sub(rownames(df.summary.fin), pattern = "'s", replacement = "")
  rownames(df.summary.fin.all) <- sub(rownames(df.summary.fin.all), pattern = "'s", replacement = "")
  if(!all(rownames(df.summary.fin.all) %in% rownames(df.summary.fin))) {stop("Creation of the clinical data summary table has failed.")}
  df.summary.fin <- merge(x = df.summary.fin, y = df.summary.fin.all, all.x = T, by.x = 0, by.y = 0, sort = F)
  df.summary.fin <- df.summary.fin[order(as.numeric(df.summary.fin[["Row.names"]] %in% sub.vars), decreasing = T),]
  rownames(df.summary.fin) <- df.summary.fin[["Row.names"]]
  df.summary.fin[["Row.names"]] <- NULL
  df.summary.fin[sub.vars, "."] <- "All_samples"
  } else {
  df <- df %>% mutate(across(.cols = c(dep.var.cats), .fns = as.factor))
  N <- nrow(df)
  df.summary <- df %>% dplyr::select(dep.var.cats, dep.var.conts, other.factors) %>% 
    sapply(., summary) %>% unlist() %>% as.data.frame()
  df.summary.fin <- rbind(df.summary, N)
  rownames(df.summary.fin)[nrow(df.summary.fin)] <- "N"
  df.summary.fin <- df.summary.fin[c(nrow(df.summary.fin), 1:(nrow(df.summary.fin)-1)), , drop = F]
}

wb <- createWorkbook()

xlsx.add.sheet <- function(df, df.name){
  if(ncol(df) > 2**14 | nrow(df) > 2**20){
    file.name.prefix <- sub(stri_reverse(substr(stri_reverse(file.name.prefix), start = 1, stop = 202)), pattern = "^\\.*", replacement = "")
    if(df.name == "Data") {
      write.table(x = df, row.names = F, quote = F, sep = ";", file = paste(file.name.prefix, df.name, "csv", sep = "."))} else {
      write.table(x = df, row.names = T, col.names = NA, quote = F, sep = ";", file = paste(file.name.prefix, df.name, "csv", sep = "."))}
} else {
  addWorksheet(wb = wb, sheetName = df.name)
  if(df.name == "Data") {
    writeData(wb = wb, sheet = df.name, x = df, rowNames = F, na.string = "NA", keepNA = T) 
  } else if(df.name == "Clin.data_summary") {
      writeData(wb = wb, sheet = df.name, x = df, rowNames = T, colNames = F, na.string = "NA", keepNA = T)
  } else {
      writeData(wb = wb, sheet = df.name, x = df, rowNames = T, na.string = "NA", keepNA = T)}}}

xlsx.add.sheet(df = df, df.name = "Data")
xlsx.add.sheet(df = df.summary.fin, df.name = "Clin.data_summary")
xlsx.add.sheet(df = res, df.name = "All_res")
xlsx.add.sheet(df = res.sig.all, df.name = "Sig_res")
xlsx.add.sheet(df = res.failed, df.name = "Failed_res")

suppressWarnings(rm(i))

if(length(wb$sheet_names) > 0) {saveWorkbook(wb = wb, file = paste(file.name.prefix, "xlsx", sep = "."), overwrite = T)}

if(nrow(res.sig) > 0) {

models <- unite(res.sig[,1:2], col = model, sep = ", subset = ")[["model"]]
times <- seq(0,1830,30)

Bs <- c(rep(100, 5), rep(50, 10), rep(25, 15), rep(15, 20), rep(10, 25), rep(5, 30), rep(3, 35), rep(1, 40))

rr.models.list <- foreach(i1 = seq(1,length(models))) %dopar% {
tryCatch(expr = {
i1 <<- i1

risk.reg.pdf <- function(m.name) {
  
  f.auc.filter <- function(auc.score) {
  n.surv <- sort(with(df.sub, eval(parse(text = sub(deparse(risk.formula), pattern = "~ 1", replacement = "")))))
  n.events <- as.numeric(grep(n.surv, pattern = "\\+", value = T, invert = T))
  min.time <- n.events[ceiling(0.1*length(n.surv))]
  if(is.na(min.time) | min.time > times[length(times)]) {min.time <- 0}
  return(subset(auc.score, subset = times >= min.time))}
  
  pdf.name.tmp <- paste("riskRegression", i1, m.name, "pdf.tmp", sep = ".")
  if(nchar(pdf.name.tmp) > 255) {
    pdf.name.tmp <- make.names(pdf.name.tmp)
    m.name <- make.names(m.name)
    pdf.name.tmp <- sub(pdf.name.tmp, pattern = m.name, replacement = "m.name")
  }
  if(nrow(sub.df)/10 <14) {width <- 14} else {width <- nrow(sub.df)/10}
    pdf(file = pdf.name.tmp, width = width)
    
    if(class(m.rr) == "Score") {
      if(any(grepl(colnames(m.rr$ROC$plotframe), pattern = "model"))) {
      if(is.multivar != any(grepl(names(m.rr$models), pattern = "Multivariate model"))) {stop("The model type evaluation failed.")}
      if(!is.null(m.rr.bs)) {
        if(any(grepl(colnames(m.rr.bs$AUC$score), pattern = "model"))) {
        if(is.multivar) {
          m.rr.bs.auc.score.multi <- subset(m.rr.bs$AUC$score, subset = model == "Multivariate model")
          m.rr.bs.auc.score.uni <- subset(m.rr.bs$AUC$score, subset = model == "Univariate model")
          if(mod.type == "Cox") {
            m.rr.bs.auc.score.multi <- f.auc.filter(m.rr.bs.auc.score.multi)
            m.rr.bs.auc.score.uni <- f.auc.filter(m.rr.bs.auc.score.uni)
            AUC.sums <- foreach(time = m.rr.bs.auc.score.uni[["times"]], .combine = c) %do% {
              m.rr.bs.auc.score.uni[m.rr.bs.auc.score.uni[["times"]] == time,][["AUC"]] + 
                m.rr.bs.auc.score.multi[m.rr.bs.auc.score.multi[["times"]] == time,][["AUC"]]
            }
            best.time.bs <- max(m.rr.bs.auc.score.uni[["times"]][which(AUC.sums == AUC.sums[which.max(AUC.sums)])])
          }
        } else {
          m.rr.bs.auc.score.uni <- subset(m.rr.bs$AUC$score, subset = model == "Univariate model")
          if(mod.type == "Cox") {
            m.rr.bs.auc.score.uni <- f.auc.filter(m.rr.bs.auc.score.uni)
            best.time.bs <- max(m.rr.bs.auc.score.uni[m.rr.bs.auc.score.uni$AUC == max(m.rr.bs.auc.score.uni$AUC[!is.na(m.rr.bs.auc.score.uni$AUC)]),][["times"]])
          }
        }
        if(is.multivar) {
          if(max(m.rr.bs.auc.score.uni$AUC[!is.na(m.rr.bs.auc.score.uni$AUC)]) <= 0.5 | max(m.rr.bs.auc.score.multi$AUC[!is.na(m.rr.bs.auc.score.multi$AUC)]) <= 0.5) {
          m.rr.auc.score.multi <- subset(m.rr$AUC$score, subset = model == "Multivariate model")
          m.rr.auc.score.uni <- subset(m.rr$AUC$score, subset = model == "Univariate model")
          if(mod.type == "Cox") {
            m.rr.auc.score.multi <- f.auc.filter(m.rr.auc.score.multi)
            m.rr.auc.score.uni <- f.auc.filter(m.rr.auc.score.uni)
            AUC.sums <- foreach(time = m.rr.auc.score.uni[["times"]], .combine = c) %do% {
              m.rr.auc.score.uni[m.rr.auc.score.uni[["times"]] == time,][["AUC"]] + 
                m.rr.auc.score.multi[m.rr.auc.score.multi[["times"]] == time,][["AUC"]]
            }
            best.time <- max(m.rr.auc.score.uni[["times"]][which(AUC.sums == AUC.sums[which.max(AUC.sums)])])
            rm(best.time.bs)
          }
        }
        } else {
          if(max(m.rr.bs.auc.score.uni$AUC[!is.na(m.rr.bs.auc.score.uni$AUC)]) <= 0.5) {
          m.rr.auc.score.uni <- subset(m.rr$AUC$score, subset = model == "Univariate model")
          if(mod.type == "Cox") {
            m.rr.auc.score.uni <- f.auc.filter(m.rr.auc.score.uni)
            best.time <- max(m.rr.auc.score.uni[m.rr.auc.score.uni$AUC == max(m.rr.auc.score.uni$AUC[!is.na(m.rr.auc.score.uni$AUC)]),][["times"]])
            rm(best.time.bs)
          }
        }
        }
      }
      } else {
        if(is.multivar) {
          m.rr.auc.score.multi <- subset(m.rr$AUC$score, subset = model == "Multivariate model")
          m.rr.auc.score.uni <- subset(m.rr$AUC$score, subset = model == "Univariate model")
          if(mod.type == "Cox") {
            m.rr.auc.score.multi <- f.auc.filter(m.rr.auc.score.multi)
            m.rr.auc.score.uni <- f.auc.filter(m.rr.auc.score.uni)
            AUC.sums <- foreach(time = m.rr.auc.score.uni[["times"]], .combine = c) %do% {
              m.rr.auc.score.uni[m.rr.auc.score.uni[["times"]] == time,][["AUC"]] + 
              m.rr.auc.score.multi[m.rr.auc.score.multi[["times"]] == time,][["AUC"]]
            }
          best.time <- max(m.rr.auc.score.uni[["times"]][which(AUC.sums == AUC.sums[which.max(AUC.sums)])])
          }
        } else {
          m.rr.auc.score.uni <- subset(m.rr$AUC$score, subset = model == "Univariate model")
          if(mod.type == "Cox") {
            m.rr.auc.score.uni <- f.auc.filter(m.rr.auc.score.uni)
            best.time <- max(m.rr.auc.score.uni[m.rr.auc.score.uni$AUC == max(m.rr.auc.score.uni$AUC[!is.na(m.rr.auc.score.uni$AUC)]),][["times"]])
          }
        }
      }
      m.rr.roc.df <- as.data.frame(m.rr$ROC$plotframe)
      if(levels(m.rr.roc.df$model)[2] == "Null model") {stop("The null model was selected for the further analysis.")}
      m.rr.roc.df.sel <- m.rr.roc.df[m.rr.roc.df[["model"]] == levels(m.rr.roc.df$model)[2],]
      m.rr.roc.df.sel[["Youden index"]] <- m.rr.roc.df.sel$TPR - m.rr.roc.df.sel$FPR
      m.rr.roc.df.sel[["Specificity"]] <- 1 - m.rr.roc.df.sel[["FPR"]]
      if(mod.type == "Cox"){
        m.rr.roc.df.sel <- m.rr.roc.df.sel[m.rr.roc.df.sel[["times"]] == if(exists("best.time.bs")) {best.time.bs} else {best.time},]}
      m.rr.roc.df.sel <- m.rr.roc.df.sel[which.max(m.rr.roc.df.sel[["Youden index"]]),]

      if(mod.type == "Cox") {
      tryCatch(expr = {
      plotAUC(m.rr, conf.int = T)
      if(exists("best.time")) {abline(v = best.time, col = "red", lty = 2, lwd = 2)}
      title(main = paste(str_wrap(models[i1], width = width * 10), "(original model(s))", sep = "\n"), cex.main = 1)},
      error = function(e){plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = paste(str_wrap(models[i1], width = width * 10), "(original model(s))", sep = "\n"), cex.main = 1, sub="ERROR: The AUC plot was not generated due to an error.")}
      )
      if(! exists("m.rr.bs")) {stop("The m.rr.bs object is missing.")}
      
        if(class(m.rr.bs) == "Score") {
          tryCatch(expr = {
          plotAUC(m.rr.bs, conf.int = T)
          if(exists("best.time.bs")) {abline(v = best.time.bs, col = "red", lty = 2, lwd = 2)}
          title(main = paste(str_wrap(models[i1], width = width * 10), paste0("(model(s) obtained after the cross-validation, based on ", m.rr.bs$split.method$B, " bootstrap samples (drawn with replacement), each of size ", m.rr.bs$split.method$M, ")"), sep = "\n"), cex.main = 1)},
          error = function(e){plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = paste(str_wrap(models[i1], width = width * 10), "(cross-validated model(s))", sep = "\n"), cex.main = 1, sub="ERROR: The AUC plot was not generated due to an error.")}
          )
          } else if(is.null(m.rr.bs)){
          plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = width * 10), cex.main = 1, sub="ERROR: The bootstrapping step was unsuccessful.")} else {
          stop("The m.rr.bs file is of unknown class.")}
      } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(1, 5), units = "null"))))
        grid.text(paste(str_wrap(models[i1], width = width * 10), "(original model(s))", sep = "\n"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1), gp = gpar(fontsize = 12, fontface = "bold"))
        grid.table(m.rr$AUC$score)
        
        if(class(m.rr.bs) == "Score") {
          grid.newpage()
          pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(1, 5), units = "null"))))
          grid.text(paste(str_wrap(models[i1], width = width * 10), paste0("(model(s) obtained after the cross-validation,\nbased on ", m.rr.bs$split.method$B, " bootstrap samples (drawn with replacement), each of size ", m.rr.bs$split.method$M, ")"), sep = "\n"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1), gp = gpar(fontsize = 12, fontface = "bold"))
          grid.table(m.rr.bs$AUC$score)  
        } else if(is.null(m.rr.bs)){
          plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = width * 10), cex.main = 1, sub="ERROR: The bootstrapping step was unsuccessful.")} else {
            stop("The m.rr.bs file is of unknown class.")}
        }
      
      if(is.multivar) {
        tryCatch(expr = {
          plotRisk(m.rr, times = if(exists("best.time.bs")) {best.time.bs} else {best.time})
          if(mod.type == "Cox"){
            title(main = paste(str_wrap(paste(models[i1], "(original models)", sep = "\n"), width = width * 10), paste("Observation time:", if(exists("best.time.bs")) {best.time.bs} else {best.time}, "days"), sep = "\n"), cex.main = 1)} else {
              title(main = str_wrap(paste(models[i1], "(original models)", sep = "\n"), width = width * 10), cex.main = 1)}},
          error = function(e){plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(paste(models[i1], "(original models)", sep = "\n"), width = width * 10), cex.main = 1, sub="ERROR: The risk plot was not generated due to an error.")})}
      
      if(is.multivar) {
        tryCatch(expr = {
          plotRisk(m.rr.bs, times = if(exists("best.time.bs")) {best.time.bs} else {best.time})
          if(mod.type == "Cox"){
            title(main = paste(str_wrap(paste(models[i1], "(models obtained after the cross-validation)", sep = "\n"), width = width * 10), paste("Observation time:", if(exists("best.time.bs")) {best.time.bs} else {best.time}, "days"), sep = "\n"), cex.main = 1)} else {
              title(main = str_wrap(paste(models[i1], "(models obtained after the cross-validation)", sep = "\n"), width = width * 10), cex.main = 1)}},
          error = function(e){plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(paste(models[i1], "(models obtained after the cross-validation)", sep = "\n"), width = width * 10), cex.main = 1, sub="ERROR: The risk plot was not generated due to an error.")})}

      tryCatch(expr = {
      plotROC(m.rr, col = cols25(), times = if(exists("best.time.bs")) {best.time.bs} else {best.time})
      points(xy.coords(x = m.rr.roc.df.sel$FPR, y = m.rr.roc.df.sel$TPR), pch=9, cex=1.5)
      abline(h = m.rr.roc.df.sel$TPR, v = m.rr.roc.df.sel$FPR, cex = 1, lty = 3)
      text(xy.coords(m.rr.roc.df.sel$FPR + 0.07, m.rr.roc.df.sel$TPR + 0.07), label=paste0("Sensitivity = ", round(m.rr.roc.df.sel[["TPR"]],3), "\nSpecificity = ",  round(m.rr.roc.df.sel[["Specificity"]], 3)))
      if(mod.type == "Cox"){
        title(main = paste(str_wrap(paste(models[i1], "(original model(s))", sep = "\n"), width = width * 10), paste("Observation time:", if(exists("best.time.bs")) {best.time.bs} else {best.time}, "days"), sep = "\n"), cex.main = 1)} else {
          title(main = str_wrap(paste(models[i1], "(original model(s))", sep = "\n"), width = width * 10), cex.main = 1)}},
      error = function(e){plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(paste(models[i1], "(original model(s))", sep = "\n"), width = width * 10), cex.main = 1, sub="ERROR: The ROC plot was not generated due to an error.")}
      )
      if(mod.type == "Cox"){
        if(is.continuous | variate == "multivar"){
          risks <- as.data.frame(predictRisk(m, newdata = sub.df, times = if(exists("best.time.bs")) {best.time.bs} else {best.time}))
          sub.df[["calculated_risks"]] <- risks[,1]
          sub.df[["risk_category"]][sub.df$calculated_risks >= m.rr.roc.df.sel$risk] <- "high"
          sub.df[["risk_category"]][sub.df$calculated_risks < m.rr.roc.df.sel$risk] <- "low"
          single.var <- gsub(sub(models[i1], pattern = "(^.*~ )([^\\+]*)( ?\\+?.*)(,.*$)", replacement = "\\2"), pattern = " *", replacement = "")
          if(is.continuous & variate == "univar"){
            colnames(sub.df)[[length(colnames(sub.df))]] <- "risk_category"
            km.formula <- as.formula(sub(as.character(res.sig[["Formula"]][i1]), pattern = "~ .*", replacement = paste("~", colnames(sub.df)[[length(colnames(sub.df))]])))
          } else {
            km.formula <- as.formula(sub(as.character(res.sig[["Formula"]][i1]), pattern = "~ .*", replacement = "~ risk_category"))
          }
          
          tryCatch(expr = {km.plot <- npsurv(data = sub.df, formula = km.formula)
          log.rank <- survdiff(data = sub.df, formula = km.formula, rho = 0)
          log.rank.p.value <- pchisq(q = log.rank$chisq, df = length(log.rank$obs)-1, lower.tail = F)}, error = function(e){plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = width * 10), cex.main = 1, sub="ERROR: The Kaplan-Meier curves generation was unsuccessful.")})} else {
          tryCatch(expr = {km.plot <- npsurv(data = sub.df, formula = formula)
          log.rank <- survdiff(data = sub.df, formula = formula, rho = 0)
          log.rank.p.value <- pchisq(q = log.rank$chisq, df = length(log.rank$obs)-1, lower.tail = F)}, error = function(e){plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = width * 10), cex.main = 1, sub="ERROR: The Kaplan-Meier curves generation was unsuccessful.")})}
        
        if(! exists("log.rank.p.value")) {log.rank.p.value <- NA} else
        if(log.rank.p.value < 0.001) {log.rank.p.value <- formatC(log.rank.p.value, format = "e", digits = 3)} else {
          log.rank.p.value <- round(log.rank.p.value,4)}
        
        if(exists("km.plot")) {
        tryCatch(expr = {
        rms::survplot(km.plot, col=cols25(), conf="none", n.risk=TRUE, dots=TRUE, lwd = 3)
        title(main = paste(str_wrap(models[i1], width = width * 10), paste("Model type:", levels(m.rr.roc.df$model)[2]), sep = "\n"), cex.main = 1)
        text(xy.coords(x = max(km.plot$time)/2, y = 1), labels = paste("log-rank test p-value:", log.rank.p.value))},
        error = function(e){plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = width * 10), cex.main = 1, sub="ERROR: The KM plot was not generated due to an error.")}
        )}
        
        tryCatch(expr = {zph <- cox.zph(m)
        zph.plot <- ggcoxzph(zph, font.main=12, caption = str_wrap(models[i1], width = width * 10))
        print(zph.plot)}, error = function(e){plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = width * 10), cex.main = 1, sub="ERROR: The Cox proportional hazard estimation returned an error.")})
      } else {
        single.var <- gsub(sub(models[i1], pattern = "(^.*~ )([^\\+]*)( ?\\+?.*)(,.*$)", replacement = "\\2"), pattern = " *", replacement = "")
      }
      if(is.continuous) {
      if(mod.type == "Cox") {
          catvar <- "risk_category"
        } else {
          catvar <- dep.var1
        }
        
      subsetting.vars <- c(single.var, catvar, unlist(str_split(dep.vars, pattern = ":"))[str_contains(x = deparse(formula), pattern = unlist(str_split(dep.vars, pattern = ":")))])
      sub.df.sub <- sub.df[rowSums(is.na(sub.df %>% dplyr::select(subsetting.vars))) == 0,]
      sub.df.sub[[catvar]] <- as.factor(sub.df.sub[[catvar]])
    
    if(nrow(sub.df.sub) > 0) {
      if(length(unique(sub.df.sub[[catvar]])) == 2) {
      mw.res <- with(sub.df.sub, wilcox.test(sub.df.sub[[single.var]] ~ sub.df.sub[[catvar]]))
      if(!is.na(mw.res$p.value)) {
        if(mw.res$p.value < 0.001) {mw.res$p.value <- formatC(mw.res$p.value, format = "e", digits = 3)} else {
          mw.res$p.value <- round(mw.res$p.value,4)}}
      test.res.label <- paste0(mw.res$method, " p-value = ", mw.res$p.value)} else {
      test.res.label <- "Wilcoxon rank sum test: independent variable has an incorrect number of levels."
      }
      box.plot <- ggplot(sub.df.sub) +
        scale_fill_manual(values = as.character(cols25())) +
 	aes(x = .data[[catvar]], y = .data[[single.var]], fill = .data[[catvar]]) +
        geom_boxplot() +
        geom_jitter(color = "azure4", size = 1) +
        stat_summary(geom = "point", fun = mean, fill = "yellow", shape = 24) +
        labs(title = paste(models[i1], if(mod.type == "Cox") {paste("(Observation time:", if(exists("best.time.bs")) {best.time.bs} else {best.time}, "days)")}, sep = "\n"),
             x = "", y = paste(single.var, "[log10 scale]", sep = "\n"), fill = catvar) +
        scale_x_discrete() + scale_y_log10()

        box.plot <- box.plot + annotate(geom = "text", x = (length(levels(box.plot$data[[catvar]]))+1)/2,
                 y = max(box.plot$data[[single.var]]) * 3,
                 label = test.res.label) +
        annotate(geom = "text", y = max(box.plot$data[[single.var]]) * 1.5, 
                 x = (length(levels(box.plot$data[[catvar]]))+1)/2,
                 label = paste("Group", paste(names(table(box.plot$data[[catvar]])), 
                                              as.vector(table(box.plot$data[[catvar]])), sep = ": n. obs. = "), collapse = "\n")) + 
        theme(plot.title = element_textbox_simple(halign = 0.5), axis.text.x=element_blank(), axis.ticks.x=element_blank())

        tryCatch(expr = {print(box.plot)},
                 error = function(e){plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = width * 10), cex.main = 1, sub="ERROR: The box plot was not generated due to an error.")})
    
        bar.plot <- ggplot(sub.df.sub) +
          geom_bar(aes(x = reorder(.data[[main.samples]], as.numeric(.data[[catvar]])), y = .data[[single.var]], fill = .data[[catvar]]), stat = "identity", show.legend = T) +
          scale_fill_manual(values = cols25()) +
          labs(title = paste(models[i1], paste0("(N = ", nrow(sub.df.sub), ")"), sep = "\n"), x = main.samples) +
          annotate(geom = "text", x = (length(box.plot$data[[catvar]])+1)/2,
                   y = max(box.plot$data[[single.var]]) * 1.05,
                   label = test.res.label) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.title = element_textbox_simple(halign = 0.5))
        tryCatch(expr = {print(bar.plot)},
                 error = function(e){plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = width * 10), cex.main = 1, sub="ERROR: The bar plot was not generated due to an error.")})
        } else {
      plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = width * 10), cex.main = 1, sub="ERROR: The boxplot and barplot generation was impossible due to the empty dataset.")
      }
      }
      } else {
        plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = width * 10), cex.main = 1, sub="ERROR: The risk estimation returned an error.")
      }
    } else {
      plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = width * 10), cex.main = 1, sub=paste(m.rr))
    }
    dev.off()
  }

  m.name <- models[i1]
  formula <- as.formula(as.character(res.sig[["Formula"]][i1]))
  risk.formula <- as.formula(sub(as.character(res.sig[["Formula"]][i1]), pattern = "~ .*", replacement = "~ 1"))
  dep.var1 <- sub(sub(deparse(risk.formula), pattern = "(Surv\\((.*)\\))? ~.*$", replacement = "\\2"), pattern = ",.*$", replacement = "")

  if(as.character(res.sig[["Data"]][i1]) == "full_table") {subs1 <- T} else {subs1 <- as.character(res.sig[["Data"]][i1])}
 
  main.var <- sub(m.name, pattern = "(^.*)( ~ )([^, ]+)([, ].*$)", replacement = "\\3")

  selected.columns <- c(main.samples, main.var, dep.var.cats, dep.var.conts, other.factors, sub.vars)
  selected.columns <- selected.columns[!is.na(selected.columns)]
  selected.columns <- selected.columns[selected.columns != "NA"]
  selected.columns <- selected.columns[selected.columns != "ALL_SAMPLES"]
  df.sub <- df %>% dplyr::select(all_of(selected.columns))

  sub.df <- subset(df.sub, subset = eval(parse(text = subs1)) & !is.na(get(dep.var1)) &
                    rowSums(as.data.frame(lapply(strsplit(sub(as.character(formula)[3], 
                                                               pattern = "^.*~ ?", replacement = ""), split = " \\+ "), 
                                                  FUN = function(x) {is.na(df.sub[,x])}))) == 0)
  
  if(grepl(m.name, pattern = "\\+")) {is.multivar <- TRUE} else {is.multivar <- FALSE}
  if(grepl(models[i1], pattern = "Surv\\(")) {
    mod.type <- "Cox"
    m <- coxph(data = sub.df, formula = as.formula(as.character(res.sig[["Formula"]][i1])), x = T)
    if(is.multivar) {
      uni.formula <<- as.formula(sub(as.character(res.sig[["Formula"]][i1]), pattern = " \\+.*$", replacement = ""))
      m.univ <- coxph(data = sub.df, formula = uni.formula, x = T)}} else {

    mod.type <- "logistic"
    m <- glm(data = sub.df, formula = as.formula(as.character(res.sig[["Formula"]][i1])), family = "binomial")
    if(is.multivar) {
    uni.formula <<- as.formula(sub(as.character(res.sig[["Formula"]][i1]), pattern = " \\+.*$", replacement = ""))
    m.univ <- glm(data = sub.df, formula = uni.formula, family = "binomial")}}
  
  test.res <- NULL
  if(is.multivar) {
    if(mod.type == "Cox"){
      m.rr <- try(riskRegression::Score(list("Multivariate model" = m, "Univariate model" = m.univ), formula = risk.formula, conf.int=TRUE, data = sub.df, times=times, metrics="auc", plots="roc", summary='risks'))
      if(class(m.rr) == "Score") {
        m.rr.bs.list <- foreach(B = Bs) %:% foreach::when(is.null(test.res)) %do% {test.res <- tryCatch(riskRegression::Score(list("Multivariate model" = m, "Univariate model" = m.univ), formula = risk.formula, conf.int=TRUE, data = sub.df, times=times, metrics="auc", plots="roc", summary='risks', split.method="bootcv", B = B, progress.bar = NULL), error = function(e){return(NULL)})}} else {
          m.rr.bs.list <- list(NULL) 
        }
      } else {
      m.rr <- try(riskRegression::Score(list("Multivariate model" = m, "Univariate model" = m.univ), formula = risk.formula, conf.int=TRUE, data = sub.df, metrics="auc", plots="roc", summary='risks'))
      if(class(m.rr) == "Score") {
      m.rr.bs.list <- foreach(B = Bs) %:% foreach::when(is.null(test.res)) %do% {test.res <- tryCatch(riskRegression::Score(list("Multivariate model" = m, "Univariate model" = m.univ), formula = risk.formula, conf.int=TRUE, data = sub.df, metrics="auc", plots="roc", summary='risks', split.method="bootcv", B = B, progress.bar = NULL), error = function(e){return(NULL)})}} else {
      m.rr.bs.list <- list(NULL)
      }
      }
    } else {
    if(mod.type == "Cox"){
      m.rr <- try(riskRegression::Score(list("Univariate model" = m), formula = risk.formula, conf.int=TRUE, data = sub.df, times=times, metrics="auc", plots="roc", summary='risks'))
      if(class(m.rr) == "Score") {
      m.rr.bs.list <- foreach(B = Bs) %:% foreach::when(is.null(test.res)) %do% {test.res <- tryCatch(riskRegression::Score(list("Univariate model" = m), formula = risk.formula, conf.int=TRUE, data = sub.df, times=times, metrics="auc", plots="roc", summary='risks', split.method="bootcv", B = B, progress.bar = NULL), error = function(e){return(NULL)})}} else {
      m.rr.bs.list <- list(NULL) 
      }
      } else {
      m.rr <- try(riskRegression::Score(list("Univariate model" = m), formula = risk.formula, conf.int=TRUE, data = sub.df, metrics="auc", plots="roc", summary='risks'))
      if(class(m.rr) == "Score") {
      m.rr.bs.list <- foreach(B = Bs) %:% foreach::when(is.null(test.res)) %do% {test.res <- tryCatch(riskRegression::Score(list("Univariate model" = m), formula = risk.formula, conf.int=TRUE, data = sub.df, metrics="auc", plots="roc", summary='risks', split.method="bootcv", B = B, progress.bar = NULL), error = function(e){return(NULL)})}} else {
      m.rr.bs.list <- list(NULL)
      }
      }
    }
      
m.rr.list <- list("m.rr" = m.rr, "m.rr.bs.list" = m.rr.bs.list)
m.rr.list$m.rr.bs.list <- m.rr.list$m.rr.bs.list[sapply(m.rr.list$m.rr.bs.list, class) == "Score"]
m.rr.bs <- if(length(m.rr.list$m.rr.bs.list) > 0) {m.rr.list$m.rr.bs.list[[1]]} else {NULL}

risk.reg.pdf(m.name = m.name)

rm(i1, uni.formula)
},
error = function(e){
  m.name <- models[i1]
  pdf.name.tmp <- paste("riskRegression", i1, m.name, "pdf.tmp", sep = ".")
  if(nchar(pdf.name.tmp) > 255) {
    pdf.name.tmp <- make.names(pdf.name.tmp)
    m.name <- make.names(m.name)
    pdf.name.tmp <- sub(pdf.name.tmp, pattern = m.name, replacement = "m.name")
  }
  pdf(file = pdf.name.tmp, width = 14)
  plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = str_wrap(models[i1], width = 14 * 10), cex.main = 1, sub=paste("ERROR:", e))
  dev.off()},
finally = suppressWarnings(rm(i1))
)
  m.rr.list
}
names(rr.models.list) <- models

pdffiles <- list.files(pattern = "^(riskRegression\\.)([0-9]+)(.*\\.pdf\\.tmp$)")
pdffiles <- pdffiles[order(as.numeric(sub(pdffiles, pattern = "^(riskRegression\\.)([0-9]+)(.*\\.pdf\\.tmp$)", replacement = "\\2")))]
pdf_combine(input = pdffiles, output = paste(file.name.prefix, "pdf", sep = "."))
unlink(pdffiles)

} else {
warning("There are no significant results of regression analyses.")
  pdf(file = paste(file.name.prefix, "pdf", sep = "."))
  plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main = "INFO: There are no significant results of regression analyses.", cex.main = 1)
  dev.off()
}
save.image(paste(file.name.prefix, "RData", sep = "."))

sessionInfo()
proc.time()
date()

cat("All done, my Friend.\n")
