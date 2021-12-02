library(phyloseq)
library(data.table)
library(tidyr)
library(biomformat)
# Plotting and analysis
library(vegan)
library(pairwiseAdonis)


# Load Metabolite Data
chem_phy <- readRDS("/Users/helenarusso/Documents/UNESP/Doutorado/UCSD/Malpighiaceae_AllSamples/Qemistree.nosync/Files_for_running_Qemistree/Subsequent analyses/ANOVAS/NEG/Malpigh_Filter_EtOH80_NEG_cladesAFGHIJ_RA.rds")


# pull out metabolite peak data
peak_data <- as.data.frame(
  as(tax_table(chem_phy),
     "matrix") )

# keep the canopus subclass, class and superclass
tax_mat <- as(tax_table(chem_phy),
              "matrix")
tax_mat <- tax_mat[ , c("canopus_subclass",
                        "canopus_class","canopus_superclass")]

sub_tax_table <- tax_table(tax_mat)

sub_phy <-phyloseq(sample_data(chem_phy),
                   otu_table(chem_phy),
                   sub_tax_table)

# sum relative areas within a class
sub_merge <- tax_glom(sub_phy,
                          "canopus_class")

# the column 'canopus_class' identifies the class
unique_classes <- unique(as(tax_table(sub_merge),
                             "matrix")[,"canopus_class"]) 

# for each class, subset the data and run anova 
class_data <- list()
for(a_class in unique_classes) {
  
  # subset compounds by class
  class_phy <- subset_taxa(sub_merge,
                             canopus_class == a_class)
  
  # create a map identifying samples by sample types (clades)
  sample_map <- setNames(sample_data(sub_phy)[["sample_type"]],
                         sample_names(sub_phy))
  sample_types <-c("cladeA",
                   "cladeF",
                   "cladeG",
                   "cladeH",
                   "cladeI",
                   "cladeJ")
  
  # get mean relative abundance for each sample type present in the class
  mean_RAs <- list() 
  for(a_sample_type in sample_types){
    mean_RAs[[a_sample_type]] <- mean(
      otu_table(
        subset_samples(class_phy,
                       sample_type == a_sample_type)))
  }
  
  # prior to modeling, normalize data by transforming relative abundances to arcsin(sqrt(x))
  arcsin_phy <- transform_sample_counts(class_phy,
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
  # output list shows for each class, mean RA by sample type, log2 fold changes, anova + tukey results
  class_data[[a_class]] <- list(
    class = a_class,
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

# put the fold changes for all classes into a dataframe
class_fold_changes <- do.call(rbind, class_data)
classes <- row.names(class_fold_changes)

class_fold_changes <-as.data.frame(apply(class_fold_changes, 2, as.numeric))
class_fold_changes$class <- classes


# adjust p values using 
class_fold_changes$adj_p_val <- p.adjust(class_fold_changes$p_val, method = "BH")

# classify enrichment in a sample type
classify_DA <- function(a_class, results_table = class_fold_changes){
  results_table <- results_table[results_table$class == a_class,]
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
class_fold_changes$sample_type_DA <-
  vapply(class_fold_changes$class,
         classify_DA,
         FUN.VALUE = character(1))


write.table(class_fold_changes,
          "/Users/helenarusso/Downloads/Malpigh_EtOH80_NEG_class_anova_and_fold_changes.tsv", sep ='\t',
          row.names = F)


