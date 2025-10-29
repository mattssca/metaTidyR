#load pacakges
#library(skimr)
#library(janitor)
#library(caret)
#library(dplyr)
#library(tibble)
#library(DataExplorer)
#library(naniar)
#library(stringr)

#set path
#setwd("../../Desktop/GIT_REPOS/metaTidyR/")

#temporarily increase max.print to a large value
options(max.print = 1e6) 

#load raw metadata
#load(file = "../../uroscanseq_analysis/data/metadata/UROSCANSEQMetadata2025_01_16.Rdata")

#create copy of raw meta
metadata = UROSCANSEQMetadata2025_01_16

#explore dataset
#skimr::skim(metadata)
#DataExplorer::create_report(metadata)
#naniar::vis_miss(metadata)

#remove leading/trailing spaces and change 'n/a' to real NA
metadata <- metadata %>%
  mutate(across(where(is.character), ~str_trim(.x))) %>%
  mutate(across(where(is.character), ~na_if(.x, "n/a"))) %>% 
  mutate(across(where(is.character), ~na_if(.x, "")))

######################################## ODD CASES #################################################


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
  dplyr::select(Status, exc0, exc2, exc3, exc4, Flag, Comment)

extra_seq_columns = metadata %>% 
  dplyr::select(Extraherat_Allprep_DNA_RNA_Blood_Mini_kit)

hist_columns = metadata %>% 
  dplyr::select(AS_Variant_Outlier,
                OHE_AS_Variant_SQ,
                OHE_AS_Variant_Src,
                OHE_AS_Variant_NE,
                OHE_AS_Variant_MP,
                OHE_AS_Variant_Ade,
                OHE_AS_Variant_Gland,
                OHE_AS_Variant_Clear,
                OHE_PrimPath_Variant_SQ,
                OHE_PrimPath_Variant_SQmetaplasia,
                OHE_PrimPath_Variant_NE,
                OHE_PrimPath_Variant_SmallCell,
                OHE_PrimPath_Variant_Ade,
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
metadata = process_column(column = "OHE_PrimPath2022_Variant_UC_with_squamous_differentiation", new_name = "primpath_sq_diff", domain = "variant_histology", type = "boolean", boolean_map = list(true = 1, false = 0)) 
metadata = process_column(column = "OHE_PrimPath2022_Variant_Sarcomatoid", new_name = "primpath_sarcomatoid", domain = "variant_histology", type = "boolean", boolean_map = list(true = 1, false = 0)) 
metadata = process_column(column = "OHE_PrimPath2022_Variant_Neuroendocrine", new_name = "primpath_neuroendocrine", domain = "variant_histology", type = "boolean", boolean_map = list(true = 1, false = 0)) 
metadata = process_column(column = "OHE_PrimPath2022_Variant_Micropapillary", new_name = "primpath_micropapilary", domain = "variant_histology", type = "boolean", boolean_map = list(true = 1, false = 0)) 
metadata = process_column(column = "OHE_PrimPath2022_Variant_UC_with_glandular_differentiation", new_name = "primpath_galnd_diff", domain = "variant_histology", type = "boolean", boolean_map = list(true = 1, false = 0)) 
metadata = process_column(column = "OHE_PrimPath2022_Variant_Clear_cell", new_name = "primpath_clear_cell", domain = "variant_histology", type = "boolean", boolean_map = list(true = 1, false = 0)) 
metadata = process_column(column = "OHE_PrimPath2022_Variant_Plasmacytoid", new_name = "primpath_plasmacytoid", domain = "variant_histology", type = "boolean", boolean_map = list(true = 1, false = 0)) 
metadata = process_column(column = "OHE_PrimPath2022_Variant_Giant_cell", new_name = "primpath_giant_cell", domain = "variant_histology", type = "boolean", boolean_map = list(true = 1, false = 0)) 
metadata = process_column(column = "OHE_PrimPath2022_Variant_Nested", new_name = "primpath_nested", domain = "variant_histology", type = "boolean", boolean_map = list(true = 1, false = 0)) 
metadata = process_column(column = "OHE_PrimPath2022_Variant_Tubular", new_name = "primpath_tubular", domain = "variant_histology", type = "boolean", boolean_map = list(true = 1, false = 0)) 
metadata = process_column(column = "Sample_date", new_name = "sample_date", domain = "clinical_data", type = "date") 
metadata = process_column(column = "Age", new_name = "age", domain = "clinical_data", type = "numeric") 
metadata = process_column(column = "EAU_Risk_Use", new_name = "is_in_set_eau_risk", domain = "cohort_info", type = "boolean", boolean_map = list(true = "Yes", false = "")) %>% mutate(is_in_set_eau_risk = tidyr::replace_na(is_in_set_eau_risk, FALSE))
metadata = process_column(column = "EAU_Risk_Over70", new_name = "eau_risk_is_over_70", domain = "clinical_data", type = "boolean", boolean_map = list(true = "Yes", false = "No"))
metadata = process_column(column = "EAU_Risk_Tumor_Status", new_name = "eau_risk_tumor_status", domain = "clinical_data", type = "factor", rename_factors = c("Primary" = "primary", "Recurrence" = "recurrence"), factor_order = c("primary", "recurrence"))
metadata = process_column(column = "EAU_Risk_Number_of_Tumors", new_name = "eau_risk_n_tumors", domain = "clinical_data", type = "factor", rename_factors = c("Multiple" = "multiple", "Single" = "single"), factor_order = c("single", "multiple"))
metadata = process_column(column = "EAU_Risk_Maximum_Tumor_Diameter", new_name = "eau_risk_tumor_size", domain = "clinical_data", type = "factor", rename_factors = c("below_3cm" = "below_3cm", "above_3cm" = "above_3cm"), factor_order = c("below_3cm", "above_3cm"))
metadata = process_column(column = "EAU_Risk_Stage", new_name = "eau_risk_tumor_stage", domain = "clinical_data", type = "factor", rename_factors = c("Ta" = "Ta", "T1" = "T1"), factor_order = c("Ta", "T1"))
metadata = process_column(column = "EAU_Risk_Concomitant_CIS", new_name = "eau_risk_is_cis", domain = "clinical_data", type = "boolean", boolean_map = list(true = "Yes", false = "No"))
metadata = process_column(column = "EAU_Risk_WHO_Grade_1973", new_name = "eau_risk_grade_who_1973", domain = "clinical_data", type = "factor", rename_factors = c("1" = "G1", "2" = "G2", "3" = "G3"), factor_order = c("G1", "G2", "G3"))
metadata = process_column(column = "EAU_Risk_VariantHist_MP_Plasma_Src_NE", new_name = "eau_risk_variant_hist", domain = "clinical_data", type = "factor", rename_factors = c("Gland" = "gland", "MP" = "mp", "NE" = "ne", "Plasma" = "plasma"), factor_order = c("gland", "mp", "ne", "plasma"))

