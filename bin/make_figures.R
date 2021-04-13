#Load packages

library(tidyverse)
library(ggExtra)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(fuzzyjoin)
library(rstudioapi)

# Set working directory to correct location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")



#' Calculate CV's within groups from Skyline summary
#' Input is name and location of .csv from document grid containing following fields:
#' @describeIn Replicate Name, Analyzer, Precursor Mz, Total Area, Peptide
#' 
#' @return CV_df, dataframe containing CVs as well as info for plotting later. Columns are:
#' Peptide: peptide sequence
#' Analyzer: mass analyzer, based on input column
#' raw_CV: % CV based on raw areas
#' median_CV: % CV based on median normalized areas
#' RT: average retention time across all reps of analyzer
#' outlier: T/F, indicating if median normalized %CV is above 20%

calc_cvs <- function(replicate_injection_results){
  
  replicate_injection_results <- read.csv(replicate_injection_results) 
  
  CV_df <- replicate_injection_results %>% group_by(Replicate.Name) %>%
    mutate(norm_area = Total.Area/median(Total.Area)) %>% ungroup() %>%
    group_by(Peptide, Analyzer) %>% 
    summarise(raw_CV = sd(Total.Area)/mean(Total.Area)*100,
              median_CV = sd(norm_area)/mean(norm_area)*100,
              RT = mean(Peptide.Retention.Time), 
              avg_area = mean(Total.Area)) %>%
    ungroup() %>%
    group_by(Analyzer) %>%
    mutate(outlier = median_CV > 20) %>% 
    ungroup() %>%
    mutate(Analyzer = gsub("m/z", "Th", Analyzer)) %>%
    mutate(Analyzer = gsub("IT", "LIT", Analyzer))
  
  return(CV_df)
  
}



#' Generate violin plot of median normalized CVs
#' CV_df from calc_cvs function
#' 
#' @return ggplot object
#' Peptides w/ %CV > 20% are plotted in gray

cv_violin <- function(CV_df) {
  
  plot <- ggplot(CV_df, aes(y = median_CV, fill = Analyzer, x = Analyzer)) + 
    geom_violin(color = NA) +
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA, color = "grey") + 
    geom_jitter(data = function(x) dplyr::filter_(x, ~ outlier), color = "grey", width = 0.1) +
    theme_minimal() + 
    scale_fill_manual(values = c("#4D908E", "#277DA1", "#F9844A", "#F9C74F")) +
    ylab("% CV") + xlab("Analyzer") +
    geom_hline(yintercept = 20, linetype = "dashed", color = "black", size = 0.5) +
    theme(legend.position = "none")
  
  return(plot)
    
}



#' Compare CV from two analyzers
#' CV_df from calc_cvs function, and the exact names of two analyzers that are to be compared
#' 
#' @return dataframe CV_pairwise, with following columns:
#' Peptide: peptide sequence
#' ratio: ratio of median norm CV in analyzer1/ median norm CV in analyzer2
#' lower_cv: Name of analyzer with lower %CV
#' log_rat: log to of ratio

cv_pairwise <- function(CV_df, analyzer1, analyzer2) {
  
  CV_pairwise <- select(CV_df, Peptide, Analyzer, median_CV) %>%
    filter(Analyzer == analyzer1 | Analyzer == analyzer2) %>%
    spread(Analyzer, median_CV) %>%
    set_colnames(c("Peptide", "a1", "a2")) %>%
    mutate(ratio = a1/a2) %>%
    mutate(lower_cv = if_else(ratio < 1, analyzer1, "Tie")) %>%
    mutate(lower_cv = if_else(ratio > 1, analyzer2, lower_cv)) %>%
    mutate(lower_cv = if_else(is.na(ratio), "Tie", lower_cv))  %>%
    mutate(log_rat = log(ratio, 10)) %>%
    select(-a1, -a2)
  
  return(CV_pairwise)
  
}



#' Plot histogram of CV ratios between two replicates calculated in cv_pairwise
#' 
#' @return Histogram of median ratios with vertical line at median ratio

cv_pairwise_histogram <- function(CV_pairwise){
  
  plot <- ggplot(CV_pairwise, aes(x = log_rat, fill = lower_cv)) +
    geom_histogram(color = "grey", bins = 50, size = 0.25) +
    geom_vline(xintercept = median(CV_pairwise$log_rat), linetype = "dashed", size = 0.5) +
    theme_minimal() +
    xlab("log10(CV in linear ion trap/CV in Orbitrap)") +
    ylab("Number of peptides") +
    scale_fill_manual(values = c("#277DA1", "#F9C74F")) +
    theme(legend.position = "bottom") +
    labs(fill = "Lower % CV in:") +
    annotate("text", label = paste0("log10(median) = ", round(median(CV_pairwise$log_rat), 3),
                                    "\nmedian = ", round(median(CV_pairwise$ratio), 3)), 
               x = median(CV_pairwise$log_rat) + 0.05,
               y = 30, hjust = "left")
  
  return(plot)
  
}



#' Compare LOQ from two analyzers
#' LIT_loq and OT_loq are outputs from calculate_loq.py 
#' 
#' @return dataframe LOQ_comp, with following columns:
#' IT_LOQ: loq in ion trap
#' Peptide: peptide sequence
#' OT_LOQ: loq in Orbitrap
#' ratio: ratio of LOQ in ion trap/ LOQ in Orbitrap
#' Lower_loq: Name of analyzer with lower LOQ
#' loq_rat: log to of ratio

comp_loq <- function(LIT_loq, OT_loq){
  LIT_fom <- read.delim(LIT_loq)
  OT_fom <- read.delim(OT_loq)
  
  LOQ_comp <- full_join(select(LIT_fom, LOQ, Peptide.Modified.Sequence), 
                      select(OT_fom, LOQ, Peptide.Modified.Sequence), by = "Peptide.Modified.Sequence") %>%
    set_colnames(c("IT_LOQ", "Peptide", "OT_LOQ")) %>% 
    mutate(ratio = IT_LOQ/OT_LOQ) %>%
    mutate(Lower_loq = if_else(ratio > 1, "Orbitrap", "Tie")) %>%
    mutate(Lower_loq = if_else(ratio < 1, "Linear ion trap", Lower_loq)) %>%
    mutate(Lower_loq = if_else(is.na(ratio), "Tie", Lower_loq)) %>%
    mutate(loq_rat = log(ratio, 10)) %>%
    filter(Lower_loq != "Tie")

  return(LOQ_comp)
}



#' Plot number of LOQs better in each analyzer
#' 
#' @return Histogram with bins for LOQ better in each analyzer

loq_count_histogram <- function(LOQ_comp){
  
  plot <- ggplot(LOQ_comp, aes(x = Lower_loq, fill = Lower_loq)) +
    geom_histogram(stat = "count") +
    theme_minimal() +
    scale_fill_manual(values = c("#277DA1", "#F9C74F")) +
    xlab("Lower LOQ")+
    coord_cartesian(ylim =c(0, 250))+
    labs(fill = "Lower LOQ in:")  +
    ylab("Number of peptides")
  
  return(plot)
  
}



#' Plot log ratio of LOQs better in each analyzer
#' 
#' @return Histogram with ratio of LOQs in each analyzer

loq_ratio_histogram <- function(LOQ_comp){
  
  plot <- ggplot(LOQ_comp, aes(x = loq_rat, fill = Lower_loq)) +
    geom_histogram(color = "grey", bins = 50, size = 0.25) +
    geom_vline(xintercept=log(median(LOQ_comp$ratio), 10), linetype = "dashed") +
    theme_minimal() +
    xlab("log10(LOQ in linear ion trap/LOQ in Orbitrap)") +
    ylab("Number of peptides") +
    scale_fill_manual(values = c("#277DA1", "#F9C74F")) +
    labs(fill = "Lower LOQ in:")  +
    annotate("text", label = paste0("log10(median) = ", round(median(LOQ_comp$loq_rat), 3),
                                    "\nmedian = ", round(median(LOQ_comp$ratio), 3)), 
               x = median(LOQ_comp$loq_rat) - 0.05,
               y = 45, hjust = "right") +
    coord_cartesian(ylim = c(0, 50))
  
  return(plot)
  
}


#' Compare LOQ before and after optimization
#' original_loq and opt_loq are outputs from calculate_loq.py
#' Should have the columns:
#' Peptide.Modified.Sequence, Protein.Name, Precursor.Mz, 
#' Precursor.Charge, LOD, LOQ, Num.Transitions
#' 
#' @return dataframe LOQ_comp
#' Peptide.Modified.Sequence: peptide sequence
#' LOQ: Calculated LOQ
#' Num.Transitions: Number of transitions for peptide
#' Transitions: Whether transitions are original or optimized

comp_loq <- function(original_loq, opt_loq){
  original_fom <- read.delim(original_loq)
  opt_fom <- read.delim(opt_loq)
  
  LOQ_comp <- mutate(select(opt_fom, Peptide.Modified.Sequence, LOQ, Num.Transitions), Transitions = "Optimized") %>%
    bind_rows(mutate(select(original_fom, Peptide.Modified.Sequence, LOQ, Num.Transitions), Transitions = "Original"))
  
  return(LOQ_comp)
}



#' Plot kernel density estimate of LOQ
#' Input results from comp_loq
#' 
#' @return ggplot of density estimate

plot_refined_loq <- function(LOQ_comp){
  original_med <- median(filter(LOQ_comp, Transitions == "Original")$LOQ)
  opt_med <- median(filter(LOQ_comp, Transitions == "Optimized")$LOQ)
  plot <-
    ggplot(LOQ_comp, aes(x = LOQ, fill = Transitions, color = Transitions))  +
    geom_density(alpha = 0, size = 1) +
    theme_minimal() +
    xlab("LOQ") +
    ylab("Number of peptides") +
    scale_fill_manual(values = c("#F94144", "#577590")) +
    scale_color_manual(values = c("#F94144", "#577590")) +
    theme(legend.position = "bottom") +
    geom_vline(xintercept = opt_med, linetype = "dashed", size = 0.5, color = "#F94144") +
    geom_vline(xintercept = original_med, linetype = "dashed", size = 0.5, color = "#577590") +
    scale_x_log10()
  
  return(plot)
}





###### Make plots



###### Plots from Figure 2

#Get AUC values from Skyline
CV_df <- calc_cvs("data/20201028_plasma_precision_info.csv")

# Violin plot
cv_violin(CV_df)

#Compare CV in ion trap to Orbitrap
CV_pairwise <- cv_pairwise(CV_df, "LIT, 1.7 Th", "OT, 1.7 Th")

#Pairwise plot
cv_pairwise_histogram(CV_pairwise)



###### Plots from Figure 3

# Compare original LOQs
original_loq <- comp_loq("data/orig_quant_limits_IT.txt", "data/orig_quant_limits_OT.txt")

# Plot histogram of number of peptides better in each histogram in original
loq_count_histogram(original_loq)

# Plot histogram of ratios of original LOQ
loq_ratio_histogram(original_loq)

# Compare optimized LOQs
opt_loq <- comp_loq("data/opt_quant_limits_IT4.txt", "data/opt_quant_limits_OT4.txt")

# Plot histogram of number of peptides better in each histogram in original
loq_count_histogram(opt_loq)

# Plot histogram of ratios of original LOQ
loq_ratio_histogram(opt_loq)


###### Plot for figure 4

# Read in data

IT_comp_LOQ <- comp_loq("data/orig_quant_limits_IT.txt", "data/opt_quant_limits_IT4.txt")

OT_comp_LOQ <- comp_loq("data/orig_quant_limits_OT.txt", "data/opt_quant_limits_OT4.txt")

# Plot for linear ion trap 
plot_refined_loq(IT_comp_LOQ)

# Plot for Orbitrap
plot_refined_loq(OT_comp_LOQ)
