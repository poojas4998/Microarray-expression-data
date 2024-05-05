library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
# library(purrr)


# ----------------------- Helper Functions to Implement ------------------------

#' Read the expression data "csv" file.
#'
#' Function to read microarray expression data stored in a csv file. The
#' function should return a sample x gene tibble, with an extra column named
#' "subject_id" that contains the geo accession ids for each subject.
#'
#' @param filename (str): the file to read.
#' @return
#' @export
#'
#' @examples expr_mat <- read_expression_table('example_intensity_data.csv')
read_expression_table <- function(filename) {
  example_intensity_data <- read_delim(filename)%>%
    pivot_longer(cols = -probe) %>%
    pivot_wider(names_from = probe) %>%
    rename("subject_id" = 1) %>%
    as_tibble()
  # options(readr.show_col_types = FALSE)
  return (example_intensity_data)
}
  #pivot_longer(, cols = -category) %>%
  # pivot_wider(, names_from = category) %>%
  #rename("category" = 1) %>%
  #as.data.frame()

#' Replaces all '.' in a string with '_'
#' @param str String to operate upon.

period_to_underscore <- function(str){
  str_output <- str_replace_all(str,"[.]", "_")
  return(str_output)
}


#' @return reformatted string.
#' @export
#'
#' @examples
#' period_to_underscore("foo.bar")
#' "foo_bar"
# period_to_underscore <- function(str) {
#  return ("")
# }


# rename variables:
# Age_at_diagnosis to Age
# SixSubtypesClassification to Subtype
# normalizationcombatbatch to Batch

#' Rename and select specified columns.
#'
#' Function to rename Age_at_diagnosis, SixSubtypesClassification, and
#' normalizationcombatbatch columns to Age, Subtype, and Batch, respectively. A
#' subset of the data should be returned, containing only the Sex, Age, TNM_Stage,
#' Tumor_Location, geo_accession, KRAS_Mutation, Subtype, and Batch columns.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) renamed and subsetted tibble
#' @export
#'
#' @examples rename_and_select(metadata)
#' 
#' 
rename_and_select <- function(data) {
  data_new <- data %>%
    rename(
           Age = Age_at_diagnosis, 
           Subtype = SixSubtypesClassification, 
           Batch = normalizationcombatbatch) %>%
    select(Sex, Age, TNM_Stage, Tumor_Location, geo_accession, KRAS_Mutation, Subtype, Batch)
  return (data_new)
}
  



#' Create new "Stage" column containing "stage " prefix.
#'
#' Creates a new column "Stage" with elements following a "stage x" format, where
#' `x` is the cancer stage data held in the existing TNM_Stage column. Stage
#' should have a factor data type.
#'
#' @param data  (tibble) metadata information for each sample
#'
#' @return (tibble) updated metadata with "Stage" column
#' @export
#'
#' @examples metadata <- stage_as_factor(metadata)
#' 
stage_as_factor <- function(df) {
  df %>% mutate(Stage = factor(paste0("stage ", TNM_Stage)))}


   #' Calculate age of samples from a specified sex.
#'
#' @param data (tibble) metadata information for each sample
#' @param sex (str) which sex to calculate mean age. Possible values are "M"
#' and "F"
#'
#' @return (float) mean age of specified samples
#' @export
#'
#' @examples mean_age_by_sex(metadata, "F")
mean_age_by_sex <- function(data, sex) {
  mean_age <- data %>% 
    filter(Sex == sex) %>% 
    summarize(mean(Age))
  return(mean_age)

}

#' Calculate average age of samples within each cancer stage. Stages should be
#' from the newly created "Stage" column.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) summarized tibble containing average age for all samples from
#' each stage.
#' @export
#'
#' @examples age_by_stage(data)
age_by_stage <- function(data) {
    avg_age_by_stage <- group_by(data,Stage)%>%
      summarise(rounded = round(mean(Age),1))
    return (avg_age_by_stage)
  }

#' Create a cross tabulated table for Subtype and Stage using dplyr methods.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) table where rows are the cancer stage of each sample, and the
#' columns are each cancer subtype. Elements represent the number of samples from
#' the corresponding stage and subtype. If no instances of a specific pair are
#' observed, a zero entry is expected.
#' @export
#'
#' @examples cross_tab <- dplyr_cross_tab(metadata)
subtype_stage_cross_tab <- function(data) {
  subtype_stage_crosstab <- data %>% dplyr::group_by(Stage) %>% count(Subtype) %>% select(-Stage) %>%
    group_by(Subtype) %>%
    pivot_wider(names_from = Subtype, values_from = n) %>% replace(is.na(.),0)
  return(subtype_stage_crosstab)
}

#' Summarize average expression and probe variability over expression matrix.
#'
#' @param exprs An (n x p) expression matrix, where n is the number of samples,
#' and p is the number of probes.
#'
#' @return A summarized tibble containing `main_exp`, `variance`, and `probe`
#' columns documenting average expression, probe variability, and probe ids,
#' respectively.
summarize_expression <- function(exprs) {
  exprs_num <- exprs %>% select(-subject_id)
  exp_summary <- gather(exprs_num, key = "probe", factor_key=TRUE)%>% group_by(probe)%>%
    summarise(mean_exp = mean(value), variance = var(value))

}

