## Load libraries --------------------------
# General utility

# TODO: get N samples in each sample type for each metabolite subclass
# TODO: come up with criteria for choosing subclasses (fold change, significance, etc.)

library(phyloseq)
library(data.table)
library(tidyr)
library(biomformat)
# Plotting and analysis
library(vegan)
library(pairwiseAdonis)

# source helper functions for plotting and subsetting
#source("src/helper_functions.R")


# set color scheme for sample types
#sample_type_cols <- c(CCA   = "#6469ed",
                     # Coral = "#e49c4c",
                      #Limu  = "#7cc854" )

# Load Cleaned Metabolite Data----------------------------------------------------------------
chem_phy         <- readRDS("/Users/helenarusso/Documents/UNESP/Doutorado/UCSD/Malpighiaceae_AllSamples/Qemistree.nosync/Files_for_running_Qemistree/Subsequent analyses/ANOVAS/NEG/Malpigh_Filter_EtOH80_NEG_cladesAFGHIJ_RA.rds")


# 13. Statistics on Metabolite Subclass by Sample Type -----------------------------------------------
# For each metabolite subclass, run anova and look for fold changes in abundance between sample types
# This information can be used to confirm that subclasss actually associate with a given sample type.
# Results inform subclass ordinations, cytoscape outputs, etc. with simple linear statistics. 

# pull out metabolite peak data
peak_data <- as.data.frame(
  as(tax_table(chem_phy),
     "matrix") )

# sum relative abundance of features in each subclass
# phyloseq treats taxonomy as hierarchical, so we need to drop most columns
# we only need 'Qemistree_sbuclass', which is the subclass
# we can also keep higher level classifications 'class' and 'superclass'
tax_mat <- as(tax_table(chem_phy),
              "matrix")
tax_mat <- tax_mat[ , c("canopus_subclass",
                        "canopus_class","canopus_superclass")]

sub_tax_table <- tax_table(tax_mat)

sub_phy <-phyloseq(sample_data(chem_phy),
                   otu_table(chem_phy),
                   sub_tax_table)

# sum relative areas within a subclass
sub_merge <- tax_glom(sub_phy,
                          "canopus_subclass")

# the column 'Qemistree_subclass' identifies the subclass
unique_subclasses <- unique(as(tax_table(sub_merge),
                             "matrix")[,"canopus_subclass"]) 

# for each subclass, subset the data and run anova 
subclass_data <- list()
for(a_subclass in unique_subclasses) {
  
  # subset compounds by subclass
  subclass_phy <- subset_taxa(sub_merge,
                             canopus_subclass == a_subclass)
  
  # create a map identifying samples by sample types
  sample_map <- setNames(sample_data(sub_phy)[["sample_type"]],
                         sample_names(sub_phy))
  sample_types <-c("cladeA",
                   "cladeF",
                   "cladeG",
                   "cladeH",
                   "cladeI",
                   "cladeJ")
  
  # get mean relative abundance for each sample type present in the subclass
  mean_RAs <- list() 
  for(a_sample_type in sample_types){
    mean_RAs[[a_sample_type]] <- mean(
      otu_table(
        subset_samples(subclass_phy,
                       sample_type == a_sample_type)))
  }
  
  # prior to modeling, normalize data by transforming relative abundances to arcsin(sqrt(x))
  arcsin_phy <- transform_sample_counts(subclass_phy,
                                        fun = function(x)
                                          asin(sqrt(x)) )
  
  # pull out transformed counts of sample_types, clean up
  sample_abunds <- data.frame(t( as(otu_table(arcsin_phy),
                                    "matrix")))
  sample_abunds[is.na(sample_abunds)] <- 0
  
  # make table with columns for sample type and abundance
  sample_table <- data.frame(sample_type = sample_map[row.names(sample_abunds)],
                             abundance = sample_abunds[[1]]
  )
 
  # count how many samples the family was observed in
  n_samp_obs <- nrow(sample_table[sample_table$abundance > 0,])
  
  # run anova on sample type and get relevant outputs
  aov_out <-  aov(abundance ~ sample_type, data = sample_table)
  aov_sum <- summary(aov_out)
  
  p_val  <- aov_sum[[1]][["Pr(>F)"]][[1]]
  sum_sq <- aov_sum[[1]]["Sum Sq"][[1,1]]
  f_val  <- aov_sum[[1]]["F value"][1,1]
  
  # run tukey HSD to test anova results
  tuk_out <-TukeyHSD(aov_out)
  
  # assemble the output as a list
  # output list shows for each subclass, mean RA by sample type, log2 fold changes, anova + tukey results
  subclass_data[[a_subclass]] <- list(
    subclass = a_subclass,
    n_samp_obs = n_samp_obs,
    cladeA_MA = mean_RAs$cladeA,
    cladeF_MA = mean_RAs$cladeF,
    cladeG_MA = mean_RAs$cladeG,
    cladeH_MA = mean_RAs$cladeH,
    cladeI_MA = mean_RAs$cladeI,
    cladeJ_MA = mean_RAs$cladeJ,
    FC_cladeF_V_cladeA = log2(mean_RAs$cladeF/mean_RAs$cladeA),
    FC_cladeG_V_cladeA = log2(mean_RAs$cladeG/mean_RAs$cladeA),
    FC_cladeH_V_cladeA = log2(mean_RAs$cladeH/mean_RAs$cladeA),
    FC_cladeI_V_cladeA = log2(mean_RAs$cladeI/mean_RAs$cladeA),
    FC_cladeJ_V_cladeA = log2(mean_RAs$cladeJ/mean_RAs$cladeA),
    FC_cladeG_V_cladeF = log2(mean_RAs$cladeG/mean_RAs$cladeF),
    FC_cladeH_V_cladeF = log2(mean_RAs$cladeH/mean_RAs$cladeF),
    FC_cladeI_V_cladeF = log2(mean_RAs$cladeI/mean_RAs$cladeF),
    FC_cladeJ_V_cladeF = log2(mean_RAs$cladeJ/mean_RAs$cladeF),
    FC_cladeH_V_cladeG = log2(mean_RAs$cladeH/mean_RAs$cladeG),
    FC_cladeI_V_cladeG = log2(mean_RAs$cladeI/mean_RAs$cladeG),
    FC_cladeJ_V_cladeG = log2(mean_RAs$cladeJ/mean_RAs$cladeG),
    FC_cladeI_V_cladeH = log2(mean_RAs$cladeI/mean_RAs$cladeH),
    FC_cladeJ_V_cladeH = log2(mean_RAs$cladeJ/mean_RAs$cladeH),
    FC_cladeJ_V_cladeI = log2(mean_RAs$cladeJ/mean_RAs$cladeI),
    tuk_cladeF_cladeA_diff  = tuk_out[[1]]["cladeF-cladeA","diff"],
    tuk_cladeF_cladeA_p     = tuk_out[[1]]["cladeF-cladeA","p adj"],
    tuk_cladeG_cladeA_diff  = tuk_out[[1]]["cladeG-cladeA","diff"],
    tuk_cladeG_cladeA_p     = tuk_out[[1]]["cladeG-cladeA","p adj"],
    tuk_cladeH_cladeA_diff  = tuk_out[[1]]["cladeH-cladeA","diff"],
    tuk_cladeH_cladeA_p     = tuk_out[[1]]["cladeH-cladeA","p adj"],
    tuk_cladeI_cladeA_diff  = tuk_out[[1]]["cladeI-cladeA","diff"],
    tuk_cladeI_cladeA_p     = tuk_out[[1]]["cladeI-cladeA","p adj"],
    tuk_cladeJ_cladeA_diff  = tuk_out[[1]]["cladeJ-cladeA","diff"],
    tuk_cladeJ_cladeA_p     = tuk_out[[1]]["cladeJ-cladeA","p adj"],
    tuk_cladeG_cladeF_diff  = tuk_out[[1]]["cladeG-cladeF","diff"],
    tuk_cladeG_cladeF_p     = tuk_out[[1]]["cladeG-cladeF","p adj"],
    tuk_cladeH_cladeF_diff  = tuk_out[[1]]["cladeH-cladeF","diff"],
    tuk_cladeH_cladeF_p     = tuk_out[[1]]["cladeH-cladeF","p adj"],
    tuk_cladeI_cladeF_diff  = tuk_out[[1]]["cladeI-cladeF","diff"],
    tuk_cladeI_cladeF_p     = tuk_out[[1]]["cladeI-cladeF","p adj"],
    tuk_cladeJ_cladeF_diff  = tuk_out[[1]]["cladeJ-cladeF","diff"],
    tuk_cladeJ_cladeF_p     = tuk_out[[1]]["cladeJ-cladeF","p adj"],
    tuk_cladeH_cladeG_diff  = tuk_out[[1]]["cladeH-cladeG","diff"],
    tuk_cladeH_cladeG_p     = tuk_out[[1]]["cladeH-cladeG","p adj"],
    tuk_cladeI_cladeG_diff  = tuk_out[[1]]["cladeI-cladeG","diff"],
    tuk_cladeI_cladeG_p     = tuk_out[[1]]["cladeI-cladeG","p adj"],
    tuk_cladeJ_cladeG_diff  = tuk_out[[1]]["cladeJ-cladeG","diff"],
    tuk_cladeJ_cladeG_p     = tuk_out[[1]]["cladeJ-cladeG","p adj"],
    tuk_cladeI_cladeH_diff  = tuk_out[[1]]["cladeI-cladeH","diff"],
    tuk_cladeI_cladeH_p     = tuk_out[[1]]["cladeI-cladeH","p adj"],
    tuk_cladeJ_cladeH_diff  = tuk_out[[1]]["cladeJ-cladeH","diff"],
    tuk_cladeJ_cladeH_p     = tuk_out[[1]]["cladeJ-cladeH","p adj"],
    tuk_cladeJ_cladeI_diff  = tuk_out[[1]]["cladeJ-cladeI","diff"],
    tuk_cladeJ_cladeI_p     = tuk_out[[1]]["cladeJ-cladeI","p adj"],
    f_stat = f_val,
    sum_of_sq = sum_sq,
    p_val = p_val
  )
  
}

# put the fold changes for all subclasss into a nice data.frame
subclass_fold_changes <- do.call(rbind, subclass_data)
subclasses <- row.names(subclass_fold_changes)

subclass_fold_changes <-as.data.frame(apply(subclass_fold_changes, 2, as.numeric))
subclass_fold_changes$subclass <- subclasses


# adjust p values using 
subclass_fold_changes$adj_p_val <- p.adjust(subclass_fold_changes$p_val, method = "BH")

# classify enrichment in a sample type
classify_DA <- function(a_subclass, results_table = subclass_fold_changes){
  results_table <- results_table[results_table$subclass == a_subclass,]
  # vector identifying enrichment in sample types 
  enriched <- c()
  # tests
  if(results_table$tuk_cladeF_cladeA_p  < 0.05 & results_table$cladeA_MA > results_table$cladeF_MA |
     results_table$tuk_cladeG_cladeA_p  < 0.05 & results_table$cladeA_MA > results_table$cladeG_MA |
     results_table$tuk_cladeH_cladeA_p  < 0.05 & results_table$cladeA_MA > results_table$cladeH_MA |
     results_table$tuk_cladeI_cladeA_p  < 0.05 & results_table$cladeA_MA > results_table$cladeI_MA |
     results_table$tuk_cladeJ_cladeA_p  < 0.05 & results_table$cladeA_MA > results_table$cladeJ_MA){
    enriched <- append(enriched,"cladeA")
  }
  if(results_table$tuk_cladeF_cladeA_p  < 0.05 & results_table$cladeF_MA > results_table$cladeA_MA |
     results_table$tuk_cladeG_cladeF_p  < 0.05 & results_table$cladeF_MA > results_table$cladeG_MA |
     results_table$tuk_cladeH_cladeF_p  < 0.05 & results_table$cladeF_MA > results_table$cladeH_MA |
     results_table$tuk_cladeI_cladeF_p  < 0.05 & results_table$cladeF_MA > results_table$cladeI_MA |
     results_table$tuk_cladeJ_cladeF_p  < 0.05 & results_table$cladeF_MA > results_table$cladeJ_MA){
    enriched <- append(enriched,"cladeF")
  }
  if(results_table$tuk_cladeG_cladeA_p  < 0.05 & results_table$cladeG_MA > results_table$cladeA_MA |
     results_table$tuk_cladeG_cladeF_p  < 0.05 & results_table$cladeG_MA > results_table$cladeF_MA |
     results_table$tuk_cladeH_cladeG_p  < 0.05 & results_table$cladeG_MA > results_table$cladeH_MA |
     results_table$tuk_cladeI_cladeG_p  < 0.05 & results_table$cladeG_MA > results_table$cladeI_MA |
     results_table$tuk_cladeJ_cladeG_p  < 0.05 & results_table$cladeG_MA > results_table$cladeJ_MA){
    enriched <- append(enriched,"cladeG")
  }
  if(results_table$tuk_cladeH_cladeA_p  < 0.05 & results_table$cladeH_MA > results_table$cladeA_MA |
     results_table$tuk_cladeH_cladeF_p  < 0.05 & results_table$cladeH_MA > results_table$cladeF_MA |
     results_table$tuk_cladeH_cladeG_p  < 0.05 & results_table$cladeH_MA > results_table$cladeG_MA |
     results_table$tuk_cladeI_cladeH_p  < 0.05 & results_table$cladeH_MA > results_table$cladeI_MA |
     results_table$tuk_cladeJ_cladeH_p  < 0.05 & results_table$cladeH_MA > results_table$cladeJ_MA){
    enriched <- append(enriched,"cladeH")
  }
  if(results_table$tuk_cladeI_cladeA_p  < 0.05 & results_table$cladeI_MA > results_table$cladeA_MA |
     results_table$tuk_cladeI_cladeA_p  < 0.05 & results_table$cladeI_MA > results_table$cladeA_MA |
     results_table$tuk_cladeI_cladeG_p  < 0.05 & results_table$cladeI_MA > results_table$cladeG_MA |
     results_table$tuk_cladeI_cladeH_p  < 0.05 & results_table$cladeI_MA > results_table$cladeH_MA |
     results_table$tuk_cladeJ_cladeI_p  < 0.05 & results_table$cladeI_MA > results_table$cladeJ_MA){
    enriched <- append(enriched,"cladeI")
  }
  if(results_table$tuk_cladeJ_cladeA_p  < 0.05 & results_table$cladeJ_MA > results_table$cladeA_MA |
     results_table$tuk_cladeJ_cladeF_p  < 0.05 & results_table$cladeJ_MA > results_table$cladeF_MA |
     results_table$tuk_cladeJ_cladeG_p  < 0.05 & results_table$cladeJ_MA > results_table$cladeG_MA |
     results_table$tuk_cladeJ_cladeH_p  < 0.05 & results_table$cladeJ_MA > results_table$cladeH_MA |
     results_table$tuk_cladeJ_cladeI_p  < 0.05 & results_table$cladeJ_MA > results_table$cladeI_MA){
    enriched <- append(enriched,"cladeJ")
  }
  if(is.null(enriched)){
    return("NA")
  }else{
    return(paste(enriched, sep = "/", collapse = ""))
  }
}

# get differentially abundant sample types for all features
subclass_fold_changes$sample_type_DA <-
  vapply(subclass_fold_changes$subclass,
         classify_DA,
         FUN.VALUE = character(1))


write.table(subclass_fold_changes,
          "/Users/helenarusso/Downloads/Malpigh_EtOH80_NEG_subclass_anova_and_fold_changes.tsv", sep ='\t',
          row.names = F)


