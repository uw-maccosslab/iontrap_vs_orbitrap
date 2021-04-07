#Load packages

library(tidyverse)
library(ggExtra)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(fuzzyjoin)
library(OrgMassSpecR)


# Set working directory to correct location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")


#' Take Skyline data about protein abundance, and generate plot of total abundance
#' Input is a list of abundances from Skyline document Grid with following fields:
#' Peptide, Total Area, Protein, Replicate, Protein Description
#' For simplicity, report was filtered to only use abundances from a single replicate
#' 
#' @return ggplot object where each point is a protein with summed abundance

protein_abundance <- function(protein_quant){
  protein_quant <- read.csv(protein_quant) %>%
    group_by(Protein, Protein.Description) %>%
    summarise(Total_area = sum(Total.Area))
  
  
  plot <- ggplot(protein_quant, aes(x=reorder(Protein, -Total_area), y = Total_area)) +
    geom_point(fill = "#264653", color = "#264653") +
    theme_minimal() + 
    theme(panel.grid.minor = element_blank(), axis.text.x = element_blank()) +
    scale_y_log10() +
    ylab("Total signal") +
    xlab("Protein") +
    coord_cartesian(ylim = c(1000000, 200000000000))
  
  return(plot)
}



#' Take Skyline info about transitions for a peptide and output plottable chromatogram
#' Input is name and location of .csv from document grid containing following fields for a single peptide:
#' @describeIn Replicate, Peptide, Fragment Ion, Raw Times, Raw Intensities, Analyzer, Min Start Time, Max End Time
#' rep1 and rep2 are Replicate names of replicates that you want to plot
#' 
#' @return chr_plot, dataframe containing intensities over time for each ion. Columns are:
#' time: time of scan (in minutes)
#' int: raw intensity of scan
#' Ion: name of ion
#' Replicate: Rep 1 or Rep 2
#' left: left integration bound
#' right: right integration bound

extract_chromatogram <- function(signal_over_time, rep1, rep2) {
  chromatogram_raw <- read.csv(signal_over_time) %>%
    filter(Replicate == rep1 | Replicate == rep2) %>%
    mutate(Raw.Times = as.character(Raw.Times)) %>%
    mutate(Raw.Intensities = as.character(Raw.Intensities)) %>%
    mutate(Replicate = if_else(Replicate == rep1, "Rep 1", "Rep 2"))
  
  chr_plot <- data.frame()

  i <- 1
  for(i in 1:nrow(chromatogram_raw)){
  info <- data.frame(str_split(chromatogram_raw$Raw.Times[i], ",")) %>%
    set_colnames(c("time")) %>% 
    cbind(data.frame(str_split(chromatogram_raw$Raw.Intensities[i], ","))) %>%
    set_colnames(c("time", "int")) %>%
    mutate(time = as.numeric(as.character(time)))  %>%
    mutate(int = as.numeric(as.character(int))) %>%
    mutate(Ion = chromatogram_raw$Fragment.Ion[i])%>%
    mutate(Replicate = chromatogram_raw$Replicate[i]) %>%
    mutate(left = chromatogram_raw$Min.Start.Time[i]) %>%
    mutate(right = chromatogram_raw$Max.End.Time[i])
  
  chr_plot <- rbind(chr_plot, info)
  }
  
return(chr_plot)
  
}



#' Generate chromatogram of all ions for peptide and reps from extract_chromatogram
#' Input is datafram returned from extract_chromatogram function
#' 
#' @return ggplot object of chromatograms
#' Vertical dashed lines show integration bounds

plot_chromatogram <- function(chr_plot){
  plot <- 
    ggplot(chr_plot, aes(x = time, y = int, color = Ion)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("#F94144", "#F9844A", "#F9C74F", "#43AA8B", "#277DA1")) +
    theme_minimal() +
    facet_wrap(~Replicate, nrow = 2) +
    geom_vline(data = chr_plot, aes(xintercept = left), linetype = "dashed") +
    geom_vline(data = chr_plot, aes(xintercept = right), linetype = "dashed") +
    ylab("Intensity") +
    xlab("time (min)") +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  
    return(plot)
}



#' Take Skyline info about transitions for a peptide and output plottable chromatogram
#' Input is name and location of .csv from document grid containing following fields for a single peptide:
#' @describeIn Replicate Name, Fragment Ion, Area
#' rep1 and rep2 are Replicate names of replicates that you want to plot
#' 
#' @return quant_df, dataframe containing info about quant for each transition and replicate. Columns are:
#' Replicate: Run 1 or Run 2
#' Ion: name of ion
#' Area: area under the curve for that ion
#' 
extract_quant <- function(quant_doc, rep1, rep2) {
 quant_df <- read.csv(quant_doc) %>%
  filter(Replicate.Name == rep1 | Replicate.Name == rep2) %>%
  select(Replicate.Name, Fragment.Ion, Area) %>%
  set_colnames(c("Replicate", "Ion", "Area")) %>%
  mutate(Replicate = if_else(Replicate == rep1, "Run 1", "Run 2"))
 
 return(quant_df)

}



#' Generate quant bar graph of each transition
#' Input is results of extract_quant function
#' 
#' @return ggplot object with bars representing area

plot_ion_quant <- function(quant_df){
  plot <- 
      ggplot(quant_df, aes(x = Replicate, fill = Ion, y = Area)) +
      geom_col()+
      scale_fill_manual(values = c("#F94144", "#F9844A", "#F9C74F", "#43AA8B", "#277DA1")) +
      theme_minimal() +
      theme(axis.title.x = element_blank(), legend.position = "right")
      
    return(plot)
}



#' Read in results of replicate injections from Skyline
#' Input is replicate_data, csv from Skyline with the following columns:
#' Replicate Name, Analyzer, Precursor Mz, Total Area, Min Start Time,
#' Max End Time, Peptide Retention Time, Total Ion Current, Peptide
#' 
#' @return replicate_info, dataframe with information about each replicate
#' Columns are:
#' Peptide: peptide sequence
#' Analyzer: description of mass analyzer
#' Precursor.Mz: mass to charge of precursor
#' total_CV: Coefficient of variance for that analyzer based on raw peak areas
#' median_CV: Median normalized CV
#' RT: average retention time
#' avg_area: average area under the curve

injection_results <- function(replicate_data){
  replicate_injection_results <- read.csv(replicate_data) 
  
  replicate_info <- replicate_injection_results %>% group_by(Replicate.Name) %>%
    mutate(norm_area = Total.Area/median(Total.Area)) %>% ungroup() %>%
    group_by(Peptide, Analyzer, Precursor.Mz) %>% 
    summarise(total_CV = sd(Total.Area)/mean(Total.Area)*100, median_CV = sd(norm_area)/mean(norm_area)*100, RT = mean(Peptide.Retention.Time), avg_area = mean(Total.Area)) %>% 
    ungroup()
  
  return(replicate_info)
}



#' Read in four dataframes from extract_metadata.py with information on each spectrum 
#' extract_metadata.py results should each have the following columns:
#' scan, precursor mz, totalIonCurrent, injectTime, scanTime
#' Also read in output from replicate_info and original metadata from Skyline (replicate_data)
#' 
#' @return results_comb, dataframe with information about each replicate. Includes:
#' Replicate.Name, Analyzer, Precursor.Mz, Total.Area, Min.Start.Time, Max.End.Time,
#' Peptide.Retention.Time, Total.Ion.Current.Area, Peptide, scan, precursormz, totalIonCurrent,
#' injectTime, scanTime, total_CV, median_CV, RT, avg_area

extract_scan_info <- function(IT_07, IT_17, OT_07, OT_17, replicate_info, replicate_data){
  
  injection_results <- read.csv(replicate_data)

  IT_07_5 <- read.csv(IT_07, row.names = 1)
  IT_07_5 <- filter(injection_results, Replicate.Name == "IT_0-7_5") %>% 
    difference_full_join(IT_07_5, by = c("Precursor.Mz" = "precursormz"), max_dist = 0.0001) %>% 
    ungroup() %>% 
    inner_join(select(filter(replicate_info, Analyzer == "IT, 0.7 m/z"), -Analyzer, -Precursor.Mz), by = "Peptide")
  
  IT_17_5 <- read.csv(IT_17, row.names = 1)
  IT_17_5 <- filter(injection_results, Replicate.Name == "IT_1-7_5") %>% 
    difference_full_join(IT_17_5, by = c("Precursor.Mz" = "precursormz"), max_dist = 0.0001) %>% 
    ungroup() %>% 
    inner_join(select(filter(replicate_info, Analyzer == "IT, 1.7 m/z"), -Analyzer, -Precursor.Mz), by = "Peptide")
  
  OT_07_5 <- read.csv(OT_07, row.names = 1)
  OT_07_5 <- filter(injection_results, Replicate.Name == "OT_0-7_5") %>% 
    difference_full_join(OT_07_5, by = c("Precursor.Mz" = "precursormz"), max_dist = 0.0001) %>% 
    ungroup() %>% 
    inner_join(select(filter(replicate_info, Analyzer == "OT, 0.7 m/z"), -Analyzer, -Precursor.Mz), by = "Peptide") %>%
    mutate(totalIonCurrent = totalIonCurrent*0.0917)
  
  OT_17_5 <- read.csv(OT_17, row.names = 1)
  OT_17_5 <- filter(injection_results, Replicate.Name == "OT_1-7_5") %>% 
    difference_full_join(OT_17_5, by = c("Precursor.Mz" = "precursormz"), max_dist = 0.0001) %>% 
    ungroup() %>% 
    inner_join(select(filter(replicate_info, Analyzer == "OT, 1.7 m/z"), -Analyzer, -Precursor.Mz), by = "Peptide") %>%
    mutate(totalIonCurrent = totalIonCurrent*0.0917)
  
  
  results_comb <- bind_rows(IT_07_5, IT_17_5, OT_07_5, OT_17_5)
  
  return(results_comb)
  
}



#' Read in results of extract_scan_info 
#' 
#' @return injection_results_summary, dataframe following columns:
#' Analyzer, Peptide, CV, total_Ions, max_Ions, mean_Ions

fill_results <- function(results_comb){
  
  injection_results_summary <- results_comb %>% 
    mutate(totalions = totalIonCurrent*injectTime/1000) %>%
    filter(scanTime >= Min.Start.Time & scanTime <= Max.End.Time) %>%
    group_by(Peptide, Replicate.Name, Analyzer) %>% 
    summarise(total_Ions = sum(totalions), mean_Ions = mean(totalions), 
                min_Ions = min(totalions), max_Ions = max(totalions), med_Ions = median(totalions),
                CV = min(median_CV))%>% 
    ungroup() %>% 
    group_by(Analyzer, Peptide) %>%
    summarise(CV = min(CV), total_Ions = mean(total_Ions), 
              max_Ions = max(max_Ions), mean_Ions = mean(mean_Ions)) %>%
    mutate(Analyzer = gsub("m/z", "Th", Analyzer)) %>%
    mutate(Analyzer = gsub("IT", "LIT", Analyzer))
  
  
  return(injection_results_summary)
  
}



#' Read in results of extract_scan_info 
#' 
#' @return IT_results, dataframe with following columns:
#' Analyzer, Peptide, CV, total_IT, max_IT, mean_IT

injecttime_results <- function(results_comb){
  
 IT_results_summary <- results_comb %>% 
    mutate(totalions = injectTime) %>%
    filter(scanTime >= Min.Start.Time & scanTime <= Max.End.Time) %>%
    group_by(Peptide, Replicate.Name, Analyzer) %>% 
   summarise(total_Ions = sum(totalions), mean_Ions = mean(totalions), 
          min_Ions = min(totalions), max_Ions = max(totalions), med_Ions = median(totalions),
          CV = min(median_CV)) %>% 
    ungroup() %>% 
    group_by(Analyzer, Peptide) %>%
    summarise(CV = min(CV), total_IT = mean(total_Ions), max_IT = max(max_Ions), 
              mean_IT = mean(mean_Ions)) %>%
    mutate(Analyzer = gsub("m/z", "Th", Analyzer)) %>%
    mutate(Analyzer = gsub("IT", "LIT", Analyzer))
  
  
  return(IT_results_summary)
  
}



#' Plot maximum fill in a single spectrum
#' Input results from  fill_results
#' 
#' @return a ggplot with the maximum ions per spectrum
#' 
plot_max_fill <- function(injection_results_summary){
  
  
  plot  <- ggplot(injection_results_summary, aes(x = Analyzer, y = max_Ions, fill = Analyzer)) +
    geom_violin(color = NA) + 
    theme_minimal() + scale_y_log10() + geom_hline(yintercept = c(10000, 4*10^5), linetype = "dashed") +
    ylab("Max ions in one scan") +
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA, color = "grey") + 
    scale_fill_manual(values = c("#4D908E", "#277DA1", "#F9844A", "#F9C74F")) +
    annotate(geom="text", x=4.4, y=12000, label="1E4",
             color="black") +
    annotate(geom="text", x=4.4, y=4.6*10^5, label="4E5",
             color="black")
  
  
  return(plot)
  
}



#' Plot mean fill in a single spectrum
#' Input results from  fill_results
#' 
#' @return a ggplot with the average ions per spectrum
#' 
plot_mean_fill <- function(injection_results_summary){
  
  
  plot  <- ggplot(injection_results_summary, aes(x = Analyzer, y = mean_Ions, fill = Analyzer)) + 
    geom_violin(color = NA) + 
    theme_minimal() + scale_y_log10() + 
    geom_hline(yintercept = c(10000, 4*10^5), linetype = "dashed") +
    ylab("Average ions in one scan") +
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA, color = "grey") + 
    scale_fill_manual(values = c("#4D908E", "#277DA1", "#F9844A", "#F9C74F")) +
    annotate(geom="text", x=4.4, y=12000, label="1E4",
             color="black") +
    annotate(geom="text", x=4.4, y=4.6*10^5, label="4E5",
             color="black")
  
  
  return(plot)
  
}



#' Plot inject time for each analyzer
#' Input results from  injecttime_results
#' 
#' @return a ggplot with the average inject time per spectrum

plot_inject_time <- function(IT_results){
  
  
  plot  <-  ggplot(IT_results, aes(x = Analyzer, y = mean_IT, fill = Analyzer)) + geom_violin(color = NA) + 
    theme_minimal() + scale_y_log10() +
    ylab("Average inject time per peptide (ms)") +
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA, color = "grey") + 
    scale_fill_manual(values = c("#4D908E", "#277DA1", "#F9844A", "#F9C74F")) 
  
  
  return(plot)
  
}



#' Predict possible fragments
#' EncyclopeDIA results, with following columns:
#' PSMId, score, q-value, posterior_error_prob, peptide, proteinIds
#' 
#' @return a dataframe of predicted fragments, all_peps_processed:
#' peptide.ms1seq: full peptide sequence
#' peptide.ms1z2: doubly charged m/z
#' peptide.ms2seq: sequence of fragment
#' peptide.ms2mz: m/z of fragment
#' peptide.ms2type: type of fragment

predict_ions <- function(EncyclopeDIA_results){
  # Read in peptides from chromatogram library and parse them
  
  encyclopedia_peps <- read.delim(EncyclopeDIA_results) %>%
    filter(q.value <= 0.01) %>%
    mutate(peptide = as.character(peptide)) %>%
    mutate(peptide_mod = gsub("-.", "", peptide)) %>%
    mutate(peptide_mod = gsub(".-", "", peptide_mod)) %>%
    mutate(peptide_mod = gsub("[[]", "", peptide_mod)) %>%
    mutate(peptide_mod = gsub("[]]", "", peptide_mod)) %>%
    mutate(peptide_mod = gsub("[+57.021464]", "", peptide_mod)) %>%
    mutate(peptide_mod = gsub("M999", "m", peptide_mod)) %>%
    select(peptide_mod)
  
  
  # Predict all fragments
  
  all_peps_processed <- apply(encyclopedia_peps, 2, FragmentPeptide, custom = list(code = "m", mass = 147.1861)) %>%
    data.frame() %>% 
    filter(grepl(']1+',peptide_mod.ms2type)) %>%
    select(peptide_mod.ms1seq, peptide_mod.ms1z2, peptide_mod.ms2seq, peptide_mod.ms2mz, peptide_mod.ms2type) %>%
    set_colnames(c("peptide.ms1seq", "peptide.ms1z2", "peptide.ms2seq", "peptide.ms2mz", "peptide.ms2type")) %>%
    filter(peptide.ms2mz >= 150 & peptide.ms2mz <= 1250)
  
  return(all_peps_processed)
}



#' Calculate number of overlapping fragments
#' Input results from predict_ions, ms1tol (width of precursor isolation window),
#' ms2tol (fragment match tolerance), tol_type(whether ms2tol is in m/z or ppm)
#' 
#' @return a dataframe of peptides and hte number of interference free transitons:
#' peptide.ms1seq: sequence
#' int_free_fragments: number of interference free fragments
#' analyzer: parameters used in simulation

ion_overlap <- function(peps_df, ms1tol, ms2tol, tol_type){
  all_ions <- peps_df
  
  
  all_ions_2 <- peps_df %>%
    set_colnames(c("pep2", "mz", "fragment", "frag_mz", "ion_type"))
  
  overlap_ms1 <- difference_full_join(all_ions, all_ions_2, by = c("peptide.ms1z2" = "mz"), max_dist = ms1tol/2) %>%
    filter(peptide.ms1seq != pep2) %>%
    mutate(delta_mass = abs(frag_mz-peptide.ms2mz)) %>%
    mutate(ppm_mass = abs((peptide.ms2mz-frag_mz)/peptide.ms2mz*10^6))
  
  if(tol_type == "delta"){
    overlap_ms1 <- overlap_ms1 %>%
      filter(delta_mass <= ms2tol) %>%
      mutate(excluded_frag = paste0(peptide.ms1seq, "_", peptide.ms2seq, "_", peptide.ms2mz))
  }
  
  if(tol_type == "ppm"){
    overlap_ms1 <- overlap_ms1 %>%
      filter(ppm_mass <= ms2tol) %>%
      mutate(excluded_frag = paste0(peptide.ms1seq, "_", peptide.ms2seq, "_", peptide.ms2mz)) %>%
      select(excluded_frag)
  }
  
  int_free_ions <- mutate(all_ions, excluded_frag = paste0(peptide.ms1seq, "_", peptide.ms2seq, "_", peptide.ms2mz)) %>%
    anti_join(overlap_ms1, by = "excluded_frag")  %>%
    filter(grepl('y',peptide.ms2type)) %>%
    select(-excluded_frag) %>%
    group_by(peptide.ms1seq) %>%
    summarise(int_free_fragments = n())
  
  return(int_free_ions)
}



#' Calculate number of overlapping fragments
#' Input results from ion_overlap for IT and OT at 1.7 Th isolation
#' and IT and OT at 0.7 m/z isolation (in that order)
#' 
#' @return a dataframe with information about interference:
#' analyzer: parameters used in simulation
#' int_free_fragments: Number of fragments for one peptide without interference
#' n: Number of peptides with that number of int_free_fragments
#' n_cumulative: Number of peptides with that number or more of int_free_fragments

int_free_total <- function(IT_1_7_sim, OT_1_7_sim, IT_0_7_sim, OT_0_7_sim){
  
  all_sim <- bind_rows(IT_1_7_sim, OT_1_7_sim, IT_0_7_sim, OT_0_7_sim)
  
  all_sim_cum <- all_sim %>%
    group_by(analyzer, int_free_fragments) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    group_by(analyzer) %>%
    arrange(-int_free_fragments) %>%
    mutate(n_cumulative = cumsum(n)) %>%
    bind_rows(data.frame(analyzer = c("LIT, 0.7 Th", "OT, 0.7 Th", "LIT, 1.7 Th", "OT, 1.7 Th"), 
                         int_free_fragments = c(0, 0, 0, 0), 
                         n = c(0, 0, 0, 0),
                         n_cumulative = c(5078, 5078, 5078, 5078))) %>%
    mutate(n_cumulative = as.integer(n_cumulative))
  
  return(all_sim_cum)
  
}



#' Plot total number of interference-free transitions
#' Input is results from int_free_total
#' 
#' @return ggplot object with number of interference free transitions

plot_total_int_free <- function(all_sim_cum){
  
  plot <- ggplot(all_sim_cum, aes(y = n_cumulative, x = int_free_fragments, color = analyzer)) +
    geom_line(size = 1) + 
    theme_minimal() +
    scale_y_log10() +
    scale_color_manual(values = c("#4D908E", "#277DA1", "#F9844A", "#F9C74F"), name = "Analyzer") +
    ylab("Number of peptides") +
    xlab("Interference-free transitions")
  
  return(plot)
  
}



#' Read in original transition set from Skyline
#' Input contains information from document grid:
#' Peptide Modified Sequence, Protein Name, Replicate Name, Precursor Mz,
#' Precursor Charge, Product Mz, Product Charge, Fragment Ion, Retention Time, 
#' Area, Total Background, Peak Rank, File Name, Isotope LabelType, Standard Type
#' 
#' @return dataframe original_trans with selected transitions:
#' fragment: Name of fragment ion
#' peptide_mod: full sequence
#' pep_ion: peptide sequence followed by _ followed by fragment name

parse_original_trans <- function(peptide_list){
  
  original_trans <- read.csv(peptide_list) %>%
    filter(Protein.Name != "Pierce standards")  %>%
    mutate(peptide_mod = gsub("[[]", "", Peptide.Modified.Sequence)) %>%
    mutate(peptide_mod = gsub("[]]", "", peptide_mod)) %>%
    mutate(peptide_mod = gsub("[+57]", "", peptide_mod)) %>%
    mutate(peptide_mod = gsub("M16", "m", peptide_mod)) %>%
    mutate(fragment = paste0("[", Fragment.Ion, "]1+")) %>%
    select(fragment, peptide_mod) %>%
    mutate(pep_ion = paste0(peptide_mod, "_", fragment))
  
  return(original_trans)
  
}



#' Calculate number of overlapping fragments
#' Input results of predict_ions, ms1tol, ms2tol, and type of tolerance
#' 
#' @return a dataframe with information on interference:
#' peptide.ms1seq, peptide.ms1z2, peptide.ms2seq, peptide.ms2mz, peptide.ms2type

ion_list <- function(peps_df, ms1tol, ms2tol, tol_type){
  all_ions <- peps_df
  
  
  all_ions_2 <- peps_df %>%
    set_colnames(c("pep2", "mz", "fragment", "frag_mz", "ion_type"))
  
  overlap_ms1 <- difference_full_join(all_ions, all_ions_2, by = c("peptide.ms1z2" = "mz"), max_dist = ms1tol/2) %>%
    filter(peptide.ms1seq != pep2) %>%
    mutate(delta_mass = abs(frag_mz-peptide.ms2mz)) %>%
    mutate(ppm_mass = abs((peptide.ms2mz-frag_mz)/peptide.ms2mz*10^6))
  
  if(tol_type == "delta"){
    overlap_ms1 <- overlap_ms1 %>%
      filter(delta_mass <= ms2tol) %>%
      mutate(excluded_frag = paste0(peptide.ms1seq, "_", peptide.ms2seq, "_", peptide.ms2mz))
  }
  
  if(tol_type == "ppm"){
    overlap_ms1 <- overlap_ms1 %>%
      filter(ppm_mass <= ms2tol) %>%
      mutate(excluded_frag = paste0(peptide.ms1seq, "_", peptide.ms2seq, "_", peptide.ms2mz)) %>%
      select(excluded_frag)
  }
  
  int_free_ions <- mutate(all_ions, excluded_frag = paste0(peptide.ms1seq, "_", peptide.ms2seq, "_", peptide.ms2mz)) %>%
    anti_join(overlap_ms1, by = "excluded_frag")  %>%
    filter(grepl('y',peptide.ms2type)) %>%
    select(-excluded_frag)
  
  return(int_free_ions)
}



#' Predict how many of the original transitions chosen were interference free
#' Output from predict_ions, output from parse_original_trans
#' 
#' @return dataframe all_sim_ions with information:
#' pep: Peptide sequence
#' n_total: total interference-free transitions for that peptide
#' n_both: transitions that were chosen and interference free for peptide
#' n_int_free: number of interference free transitions
#' chosen: number of chosen transitions
#' frac_chosen_int_free: portion of chosen transitions w/o interference
#' analyzer: analyzer description

eval_original_trans <- function(all_peps_processed, parse_original_trans){
  
  IT_1_7_sim_ions <- ion_list(all_peps_processed, 1.7, 0.7, "delta") %>%
    mutate(pep_ion = paste0(peptide.ms1seq, "_", peptide.ms2type)) %>%
    semi_join(original_trans, by= c("peptide.ms1seq" = "peptide_mod")) %>%
    full_join(original_trans, by = "pep_ion") %>%
    distinct(.keep_all = TRUE) %>%
    mutate(int_free = if_else(is.na(peptide.ms1seq), 0, 1)) %>%
    mutate(chosen = if_else(is.na(peptide_mod), 0, 1)) %>%
    mutate(pep = gsub("\\_.*", "", pep_ion)) %>%
    select(pep, int_free, chosen) %>%
    mutate(both = int_free + chosen) %>%
    mutate(both = if_else(both == 2, 1, 0)) %>%
    group_by(pep) %>%
    summarise(n_total = n(), n_both = sum(both), n_int_free = sum(int_free), chosen = sum(chosen)) %>%
    mutate(frac_chosen_int_free = n_both/chosen) %>%
    mutate(analyzer = "LIT, 1.7 Th")
  
  OT_1_7_sim_ions <- ion_list(all_peps_processed, 1.7, 20, "ppm") %>%
    mutate(pep_ion = paste0(peptide.ms1seq, "_", peptide.ms2type)) %>%
    semi_join(original_trans, by= c("peptide.ms1seq" = "peptide_mod")) %>%
    full_join(original_trans, by = "pep_ion") %>%
    distinct(.keep_all = TRUE) %>%
    mutate(int_free = if_else(is.na(peptide.ms1seq), 0, 1)) %>%
    mutate(chosen = if_else(is.na(peptide_mod), 0, 1)) %>%
    mutate(pep = gsub("\\_.*", "", pep_ion)) %>%
    select(pep, int_free, chosen) %>%
    mutate(both = int_free + chosen) %>%
    mutate(both = if_else(both == 2, 1, 0)) %>%
    group_by(pep) %>%
    summarise(n_total = n(), n_both = sum(both), n_int_free = sum(int_free), chosen = sum(chosen)) %>%
    mutate(frac_chosen_int_free = n_both/chosen) %>%
    mutate(analyzer = "OT, 1.7 Th")
  
  IT_0_7_sim_ions <- ion_list(all_peps_processed, 0.7, 0.7, "delta") %>%
    mutate(pep_ion = paste0(peptide.ms1seq, "_", peptide.ms2type)) %>%
    semi_join(original_trans, by= c("peptide.ms1seq" = "peptide_mod")) %>%
    full_join(original_trans, by = "pep_ion") %>%
    distinct(.keep_all = TRUE) %>%
    mutate(int_free = if_else(is.na(peptide.ms1seq), 0, 1)) %>%
    mutate(chosen = if_else(is.na(peptide_mod), 0, 1)) %>%
    mutate(pep = gsub("\\_.*", "", pep_ion)) %>%
    select(pep, int_free, chosen) %>%
    mutate(both = int_free + chosen) %>%
    mutate(both = if_else(both == 2, 1, 0)) %>%
    group_by(pep) %>%
    summarise(n_total = n(), n_both = sum(both), n_int_free = sum(int_free), chosen = sum(chosen)) %>%
    mutate(frac_chosen_int_free = n_both/chosen) %>%
    mutate(analyzer = "LIT, 0.7 Th")
  
  OT_0_7_sim_ions <- ion_list(all_peps_processed, 0.7, 20, "ppm") %>%
    mutate(pep_ion = paste0(peptide.ms1seq, "_", peptide.ms2type)) %>%
    semi_join(original_trans, by= c("peptide.ms1seq" = "peptide_mod")) %>%
    full_join(original_trans, by = "pep_ion") %>%
    distinct(.keep_all = TRUE) %>%
    mutate(int_free = if_else(is.na(peptide.ms1seq), 0, 1)) %>%
    mutate(chosen = if_else(is.na(peptide_mod), 0, 1)) %>%
    mutate(pep = gsub("\\_.*", "", pep_ion)) %>%
    select(pep, int_free, chosen) %>%
    mutate(both = int_free + chosen) %>%
    mutate(both = if_else(both == 2, 1, 0)) %>%
    group_by(pep) %>%
    summarise(n_total = n(), n_both = sum(both), n_int_free = sum(int_free), chosen = sum(chosen)) %>%
    mutate(frac_chosen_int_free = n_both/chosen) %>%
    mutate(analyzer = "OT, 0.7 Th")
  
  all_sim_ions <- bind_rows(IT_1_7_sim_ions, OT_1_7_sim_ions, IT_0_7_sim_ions, OT_0_7_sim_ions)
  
  return(all_sim_ions)
  
}



#' Generate violin plot of selected interference-free transitions per peptide
#' Input results from eval_original_trans
#' 
#' @return ggplot object showing how good initial transition set was

plot_eval_original_trans <- function(all_sim_ions){
  
  plot <- ggplot(all_sim_ions, aes(y = n_both, fill = analyzer, x = analyzer)) +
    geom_violin(color = NA, adjust =2) +
    theme_minimal() +
    ylab("Interference-free transitions chosen originally") +
    scale_fill_manual(values = c("#4D908E", "#277DA1", "#F9844A", "#F9C74F")) +
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA, color = "grey") +
    theme(legend.position = "none")

  
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



#' Plot histogram of number of transitions before and after refinement
#' Input results from comp_loq
#' 
#' @return ggplot of transition count

plot_refined_transitions <- function(LOQ_comp){
  plot <-
    ggplot(IT_comp_LOQ, aes(x = Num.Transitions, fill = Transitions))  +
    geom_histogram(color = "white", position = "dodge", binwidth = 1) +
    theme_minimal() +
    xlab("Number of transitions") +
    ylab("Number of peptides") +
    scale_fill_manual(values = c("#F94144", "#577590")) +
    theme(legend.position = "bottom") +
    coord_cartesian(xlim = c(3, 11)) +
    scale_x_continuous(breaks = c(4, 6, 8, 10))
  
  return(plot)
}



###### Make figures



###### Plot from Figure S2

protein_abundance("data/protein_quant_total.csv")



###### Plots from Figure S4

# Read in chromatogram info
chr_plot <- extract_chromatogram("data/interference_signal_over_time.csv", "IT_1-7_3", "IT_1-7_6")

# Plot chromatogram
plot_chromatogram(chr_plot)

#Read in quant info
quant_df <- extract_quant("data/interference_transition_quant.csv", "IT_1-7_3", "IT_1-7_6")

#Plot quant info
plot_ion_quant(quant_df)



###### Plots from Figure S5

#Parse results
replicate_info <- injection_results("data/20201028_plasma_precision_info.csv") 

results_summary <- extract_scan_info("data/quant_plasma_IT_0-7_5.csv", 
                                     "data/quant_plasma_IT_1-7_5.csv", 
                                     "data/quant_plasma_OT_0-7_5.csv",
                                     "data/quant_plasma_OT_1-7_5.csv",
                                     replicate_info,
                                     "data/20201028_plasma_precision_info.csv")

fill_results_summary <- fill_results(results_summary)

IT_results_summary <- injecttime_results(results_summary)

# Violin plot of maximum ions per spectrum
plot_max_fill(fill_results_summary)

# Violin plot of mean ions per spectrum
plot_mean_fill(fill_results_summary)

# Violin plot of average inject time
plot_inject_time(IT_results_summary)



##### Plot from Figure S6 

# Read in EncyclopeDIA results
all_peps_processed <- predict_ions("data/20191014_LRH_lumos_GW_lib_concatenated_results.txt")

#Predict interference under each condition
IT_1_7_sim <- ion_overlap(all_peps_processed, 1.7, 0.7, "delta") %>%
  mutate(analyzer = "LIT, 1.7 Th")

OT_1_7_sim <- ion_overlap(all_peps_processed, 1.7, 20, "ppm") %>%
  mutate(analyzer = "OT, 1.7 Th")

IT_0_7_sim <- ion_overlap(all_peps_processed,0.7, 0.7, "delta") %>%
  mutate(analyzer = "LIT, 0.7 Th")

OT_0_7_sim <- ion_overlap(all_peps_processed, 0.7, 20, "ppm") %>%
  mutate(analyzer = "OT, 0.7 Th")

all_sim_cum <- int_free_total(IT_1_7_sim, OT_1_7_sim, IT_0_7_sim, OT_0_7_sim)

# Generate plot 
plot_total_int_free(all_sim_cum)



##### Plot from Figure S7

#Parse original transition set
original_trans <- parse_original_trans("data/transition_quant_IT_original.csv")

# Determine which transitions were selected as interference free from original set 
all_sim_ions <- eval_original_trans(all_peps_processed, parse_original_trans)

# Make plot
plot_eval_original_trans(all_sim_ions)



###### Plot for figure S8

# Read in data

IT_comp_LOQ <- comp_loq("data/orig_quant_limits_IT.txt", "data/opt_quant_limits_IT4.txt")

OT_comp_LOQ <- comp_loq("data/orig_quant_limits_OT.txt", "data/opt_quant_limits_OT4.txt")

# Plot for linear ion trap 
plot_refined_loq(IT_comp_LOQ)

# Plot for Orbitrap
plot_refined_loq(OT_comp_LOQ)



###### Plot for figure S9

# Plot for linear ion trap 
plot_refined_transitions(IT_comp_LOQ)

# Plot for Orbitrap
plot_refined_transitions(OT_comp_LOQ)


