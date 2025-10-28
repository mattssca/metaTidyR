#load pacakges
library(skimr)
library(janitor)
library(caret)
library(dplyr)
library(tibble)
library(DataExplorer)
library(naniar)
library(stringr)

#set path
setwd("~/BIOINFORMATICS/metaTidyR/")

#temporarily increase max.print to a large value
options(max.print = 1e6) 

#load raw metadata
load(file = "data/raw/UROSCANSEQMetadata2025_01_16.Rdata")

#create copy of raw meta
metadata = UROSCANSEQMetadata2025_01_16

#explore dataset
skimr::skim(metadata)
DataExplorer::create_report(metadata)
naniar::vis_miss(metadata)

#remove leading/trailing spaces and change 'n/a' to real NA
metadata <- metadata %>%
  mutate(across(where(is.character), ~str_trim(.x))) %>%
  mutate(across(where(is.character), ~na_if(.x, "n/a"))) %>% 
  mutate(across(where(is.character), ~na_if(.x, "")))

######################################## ODD CASES #################################################
################################## AS PRIMEPATH VARIANTS ###########################################

#Select the columns to combine
#as pathology
as_columns <- metadata %>%
  select(OHE_AS_Variant_SQ:OHE_AS_Variant_Clear)

as_columns <- as_columns %>%
  mutate(as_variants = apply(., 1, function(row) {
    present_variants <- names(row)[which(row == 1)]
    if (length(present_variants) > 0) {
      tolower(paste(present_variants, collapse = ", "))
    } else {
      NA
    }
  }))

as_variants = as.data.frame(as_columns$as_variants)

#primepath patology
prime_path_columns <- metadata %>%
  select(OHE_PrimPath_Variant_SQ:OHE_PrimPath2022_Variant_Tubular)

prime_path_columns <- prime_path_columns %>%
  mutate(prime_path_variants = apply(., 1, function(row) {
    present_variants <- names(row)[which(row == 1)]
    if (length(present_variants) > 0) {
      tolower(paste(present_variants, collapse = ", "))
    } else {
      NA
    }
  }))

prime_path_variants = as.data.frame(prime_path_columns$prime_path_variants)

#add to metadata
#tmp shift rownames
metadata <- metadata %>%
  rownames_to_column(var = "rowname")

as_variants <- as_variants %>%
  rownames_to_column(var = "rowname")

prime_path_variants <- prime_path_variants %>%
  rownames_to_column(var = "rowname")

#perform the left join using the "rowname" column
metadata <- metadata %>%
  left_join(as_variants, by = "rowname") %>% 
  rename(as_variants = `as_columns$as_variants`)

metadata <- metadata %>%
  left_join(prime_path_variants, by = "rowname") %>% 
  rename(prime_path_variants = `prime_path_columns$prime_path_variants`)

#remove the "rowname" column and restore row names
metadata <- metadata %>%
  column_to_rownames(var = "rowname")

####################################### HIST COLUMN ################################################
#standardize the column
hist_meta = metadata %>% 
  select(Multiple_Hist_Subtype) %>% 
  mutate(hist_sub = str_to_lower(Multiple_Hist_Subtype),
         hist_sub = str_trim(hist_sub))
  
hist_meta <- hist_meta %>%
  mutate(hist_sub = ifelse(Multiple_Hist_Subtype == "N", "none", hist_sub))

subtypes = c("micropapillary", "nested", "large_nested", "tubular_microcystic", 
             "plasmacytoid", "sarcomatoid", "lipid_rich", "lymphoepithelioma_like", 
             "clear_cell_glycogene_rich", "giant_cell", "poorly_differentiated")

# Named vector for mapping
subtype_map <- c(
  "micropapillär" = "micropapillary",
  "mikropapillär" = "micropapillary",
  "micropapillär + jättecellsmorfologi" = "micropapillary_giant",
  "micropapillär + skivepitel" = "micropapillary",
  "mikropapillär & nested" = "micropapillary",
  "nested" = "nested",
  "nested variant" = "nested",
  "large_nested" = "large_nested",
  "tubulär och skivepitekdifferentiering" = "tubular_microcystic",
  "tubular_microcystic" = "tubular_microcystic",
  "plasmacytoid" = "plasmacytoid",
  "plasmacytoid + sarkomatoid" = "plasmacytoid",
  "sarkomatoid" = "sarcomatoid",
  "sarcomatoid" = "sarcomatoid",
  "lipid_rich" = "lipid_rich",
  "lymphoepithelioma_like" = "lymphoepithelioma_like",
  "klarcellig differentiering" = "clear_cell_glycogene_rich",
  "clear_cell_glycogene_rich" = "clear_cell_glycogene_rich",
  "giant_cell" = "giant_cell",
  "pleomorphic giant cell carcinoma" = "giant_cell",
  "poorly_differentiated" = "poorly_differentiated",
  "neuroendokrin" = "poorly_differentiated",
  "neuroendokrin differentiering" = "poorly_differentiated",
  "storcelllig" = "poorly_differentiated",
  "småcellig" = "poorly_differentiated",
  "adenodifferntiering" = "poorly_differentiated",
  "glandulär differentiering" = "poorly_differentiated",
  "glandulärdifferentiering" = "poorly_differentiated",
  "körteldifferentiering" = "poorly_differentiated",
  "mucinös och signetringcellsdifferentiering" = "poorly_differentiated",
  "signetringscellkomponent + rhabdoid komponent" = "poorly_differentiated",
  "flera inkl. småcellig variant" = "poorly_differentiated",
  "skivepitelcancer" = "poorly_differentiated",
  "skivepiteldiff" = "poorly_differentiated",
  "skivepiteldifferentiering" = "poorly_differentiated",
  "skivepitelmetaplasi" = "poorly_differentiated",
  "nej" = NA,
  "n" = NA,
  "none" = "none",
  "j" = NA,
  "m" = NA
)

hist_meta <- hist_meta %>%
  mutate(consensus_subtype = subtype_map[hist_sub])
  
###################################### REMOVE COLUMNS ##############################################
#get all columns
all_columns = colnames(metadata)

#subset columns that are to be removed
remove_columns = metadata %>% 
  dplyr::select(RNA_cohort_name,
                CategoryGroup,
                TMA_OnTMA,
                exc1,
                Seq_name, 
                RNA_sequenced,
                XRNA_cohort_name.1, 
                Labb_ID, 
                EAU_Risk_Score_Class, 
                EAU_Risk_Score_Class_Adjusted, 
                TMA_SampleID,
                UROSCANSEQ_COHORT,
                Predictions_5classes,
                Predictions_7classes,
                Uro, UroA, UroB, UroC, GU, BaSq, Mes, ScNE, 
                trialno)

comment_columns = metadata %>% 
  select(Status, exc0, exc2, exc3, exc4, Flag, Comment)

extra_seq_columns = metadata %>% 
  select(Extraherat_Allprep_DNA_RNA_Blood_Mini_kit)

hist_columns = metadata %>% 
  select(OHE_AS_Variant_SQ,
         OHE_AS_Variant_Src,
         OHE_AS_Variant_NE,
         OHE_AS_Variant_MP,
         OHE_AS_Variant_Ade,
         OHE_AS_Variant_Gland,
         OHE_AS_Variant_Clear,
         AS_Variant_Outlier,
         OHE_PrimPath_Variant_SQ,
         OHE_PrimPath_Variant_SQmetaplasia,
         OHE_PrimPath2022_Variant_UC_with_squamous_differentiation,
         OHE_PrimPath2022_Variant_Sarcomatoid,
         OHE_PrimPath_Variant_NE,
         OHE_PrimPath_Variant_SmallCell,
         OHE_PrimPath2022_Variant_Neuroendocrine,
         OHE_PrimPath2022_Variant_Micropapillary,
         OHE_PrimPath_Variant_Ade,
         OHE_PrimPath2022_Variant_UC_with_glandular_differentiation,
         OHE_PrimPath2022_Variant_Clear_cell,
         OHE_PrimPath2022_Variant_Plasmacytoid,
         OHE_PrimPath2022_Variant_Giant_cell,
         OHE_PrimPath2022_Variant_Nested,
         OHE_PrimPath2022_Variant_Tubular, 
         Multiple_Hist_Subtype)

#subtract these columns from the main metadata
metadata <- metadata[, setdiff(names(metadata), names(remove_columns))]
metadata <- metadata[, setdiff(names(metadata), names(comment_columns))]
metadata <- metadata[, setdiff(names(metadata), names(extra_seq_columns))]
metadata <- metadata[, setdiff(names(metadata), names(hist_columns))]

###################################### COLUMN BY COLUMN ############################################
#source helper function
source("R/process_columns.R")

#initialize the column change log
column_change_log <- list()

metadata = process_column(column = "XRNA_cohort_name", new_name = "sample_id", domain = "identifier")
metadata = process_column(column = "Keep_Remove", new_name = "keep", domain = "qc_fields", type = "boolean", boolean_map = list(true = "Keep", false = "Remove"))
metadata = process_column(column = "In719set", new_name = "is_in_set_719", domain = "cohort_info", type = "boolean", boolean_map = list(true = 1, false = 0))
metadata = process_column(column = "In676set", new_name = "is_in_set_676", domain = "cohort_info", type = "boolean", boolean_map = list(true = 1, false = 0))
metadata = process_column(column = "InHQ_572set", new_name = "is_in_set_hq_572", domain = "cohort_info", type = "boolean", boolean_map = list(true = 1, false = 0))
metadata = process_column(column = "InHQ_572set_IndexUC_533", new_name = "is_in_set_hq_572_index_uc_533", domain = "cohort_info", type = "boolean", boolean_map = list(true = 1, false = 0))
metadata = process_column(column = "InLQ_147set", new_name = "is_in_set_lq_147", domain = "cohort_info", type = "boolean", boolean_map = list(true = 1, false = 0))
metadata = process_column(column = "InLQ_147set_IndexUC_129", new_name = "is_in_set_lq_147_index_uc_129", domain = "cohort_info", type = "boolean", boolean_map = list(true = 1, false = 0))
metadata = process_column(column = "QC_removal", new_name = "qc_comment", domain = "qc_fields", type = "factor", factor_levels = c("Bad_batches", "High_Quality", "Low_Quality"), rename_factors = c("Bad_batches" = "bad_batch", "High_Quality" = "high_quality", "Low_Quality" = "low_quality"), factor_order = c("bad_batch", "low_quality", "high_quality")) 
metadata = process_column(column = "Category", new_name = "sample_category", domain = "cohort_info", type = "factor", rename_factors = c("Non_UC" = "non_uc", "Recurrence" = "recurrence", "Replicate" = "replicate", "ReSeq" = "reseq", "UC_index" = "uc_index")) 
metadata = process_column(column = "IsOnTMA", new_name = "is_in_set_tma", domain = "cohort_info", type = "boolean", boolean_map = list(true = "YES", false = "NO"))
metadata = process_column(column = "in_RNA_TMA_set", new_name = "is_in_set_rna_tma", domain = "cohort_info", type = "boolean", boolean_map = list(true = 1, false = 0))
metadata = process_column(column = "arm", new_name = "is_in_set_arm_tma", domain = "cohort_info", type = "boolean", boolean_map = list(true = "On TMA", false = "Not on TMA"))
metadata = process_column(column = "SequencingBatch", new_name = "seq_batch", domain = "seq_info", type = "factor", rename_factors = c("CTGbatch024" = "ctg_batch_024", "CTGbatch026" = "ctg_batch_026", "CTGbatch110" = "ctg_batch_110", "CTGbatch121" = "ctg_batch_121", "CTGbatch168" = "ctg_batch_168", "CTGbatch2022_176" = "ctg_batch_176")) 
metadata = process_column(column = "CMDvsCTG", new_name = "seq_facility", domain = "seq_info", type = "factor", rename_factors = c("CMD" = "cmd", "CTG" = "ctg")) 
metadata = process_column(column = "instrument", new_name = "seq_instrument", domain = "seq_info", type = "factor", rename_factors = c("@A00681" = "A00681", "@NB501697" = "NB501697", "@NB501699" = "NB501699")) 
metadata = process_column(column = "run_number", new_name = "run_number", domain = "seq_info", type = "numeric") 
metadata = process_column(column = "flowcell_ID", new_name = "flowcell_id", domain = "seq_info") 
metadata = process_column(column = "DNA_Nanodrop_ng_ul", new_name = "dna_nanodrop_ng_ul", domain = "qc_fields", type = "double") 
metadata = process_column(column = "DNA_260_280", new_name = "dna_260_280", domain = "qc_fields", type = "double") 
metadata = process_column(column = "DNA_260_230", new_name = "dna_260_230", domain = "qc_fields", type = "double") 
metadata = process_column(column = "RNA_Nanodrop_ng_ul", new_name = "rna_nanodrop_ng_ul", domain = "qc_fields", type = "double") 
metadata = process_column(column = "RNA_260_280", new_name = "rna_260_280", domain = "qc_fields", type = "double") 
metadata = process_column(column = "RNA_260_230", new_name = "rna_260_230", domain = "qc_fields", type = "double") 
metadata = process_column(column = "Dilution", new_name = "dilution", domain = "qc_fields") %>% mutate(dilution = recode(dilution, "Clarity" = "clarity", "Intorkad" = "dried", "ospätt" = "non_diluted", "Stor bit" = "big_sample"))
metadata = process_column(column = "RNA_HS_Qubit_ng_ul", new_name = "rna_hs_qubit_ng_ul", domain = "qc_fields", type = "double") 
metadata = process_column(column = "RIN_value", new_name = "rin_value", domain = "qc_fields", type = "double") 
