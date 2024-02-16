#' Impute missing peaks with provided feature-level imputation values
#'
#'@param mzroll_list: data in triple omic structure
#'@param lod_values: a tibble that maps groupId to log2 feature-level imputation values
#' if a tibble is not provided, half min value per feature will be used as imputation values
#'@param quant_var: column to use for peak values, the function performs imputation in log2 space, peak quant values must be in log2 space
#'@param imp_sd: standard deviation of Gaussian distribution to use for missing peak imputation
#'
#'@return triple omic data with imputed missing peaks
#'
#'@export
impute_missing_peaks <- function(mzroll_list,
                                 lod_values = NULL,
                                 quant_var = "log2_abundance",
                                 imp_sd = 0.15) {

  claman::test_mzroll_list(mzroll_list)

  ##  get half min value per feature to use for imputation if feature-specific imputation values are not provided
  if(is.null(lod_values)) {
    lod_values <- mzroll_list$measurements %>%
      group_by(groupId) %>%
      filter(!!rlang::sym(quant_var) == min(!!rlang::sym(quant_var))/2) %>%
      ungroup() %>%
      select(groupId, !!rlang::sym(quant_var))
  }

  ## check that required columns are present in the imputation tibble
  stopifnot(colnames(lod_values) %in% c("groupId", rlang::sym(quant_var)))

  ## check if imputation value is unique per feature
  if (nrow(lod_values) > nrow(lod_values %>% dplyr::distinct(groupId, .keep_all = TRUE))) {
    stop("only one value per feature must be provided to impute missing peaks")
  }

  features <- mzroll_list$features %>% dplyr::select(groupId)

  ## check if groupIds match between the data and imputation tibble
  if(!all(features$groupId %in% lod_values$groupId)) {
    stop("groupId values must match between lod_value df and feature table of triple omic data")
  }

  valid_quant_var <- setdiff(
    mzroll_list$design$measurements$variable,
    c(mzroll_list$design$feature_pk, mzroll_list$design$sample_pk)
  )

  checkmate::assertChoice(quant_var, valid_quant_var)
  checkmate::assertNumeric(mzroll_list$measurements[[quant_var]])

  ## find missing peaks
  missing_peaks <- tidyr::expand_grid(
    groupId = mzroll_list$features$groupId,
    sampleId = mzroll_list$samples$sampleId) %>%
    dplyr::anti_join (
      mzroll_list$measurements,
      by = c("groupId", "sampleId")
    )

  ## impute missing peaks
  missing_peaks_imputed <- dplyr::left_join(
    missing_peaks,
    lod_values,
    by = c("groupId")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(!!rlang::sym(quant_var) := stats::rnorm(1, mean = !!rlang::sym(quant_var)+1, sd = imp_sd))

  ## merge measured peaks with imputed peaks
  completed_peaks <- dplyr::bind_rows(
    mzroll_list$measurements,
    missing_peaks_imputed
  )

  mzroll_list$measurements <- completed_peaks
  return(mzroll_list)
}

#' Find comparisons that have at least one imputed value
#'
#'@param feature_id: groupId for specific feature
#'@param cond_num: ConditionNum in metadata
#'@param metadata: df of metadata
#'@param df: data in long format
#'
#'@export
imputed_comparisons <- function(feature_id, cond_num, metadata, df) {

  output <- data.frame(
    term = character(1),
    groupId = factor(1),
    imputed = "no")

  if (nrow(metadata %>% dplyr::filter(ConditionNum == cond_num)) < 2)
    {return(output)}
  else {
    ref_cond_num <- unique(metadata$ReferenceConditionNum[metadata$ConditionNum == cond_num])
    cond <- unique(metadata$condition[metadata$ConditionNum == cond_num])
    ref_cond <- unique(metadata$condition[metadata$ConditionNum == ref_cond_num])

    ## skip if condition and reference condition are the same
    if (cond == ref_cond) {return(output)} else {
      df <- df %>%
        dplyr::filter(ConditionNum %in% c(cond_num, ref_cond_num)) %>%
        dplyr::mutate(condition = factor(condition, levels = c(ref_cond, cond))) %>%
        dplyr::filter(groupId == feature_id) %>%
        dplyr::select(imputed)
      if ("yes" %in% df$imputed) {
        output$term = paste("condition", cond, sep = "")
        output$groupId = feature_id
        output$imputed = "yes"
        }
      }
    return(output)
    }
  }

#' Generate complete dataset upon imputing missing peaks
#'
#'@param mzroll_list: data in triple omic structure that's already filtered to only keep peakgroups of interest
#'@param metadata: metadata of experiment
#'@param lod_values: a tibble that maps groupId to log2 feature-level imputation values
#'
#'@return dataframe of complete dataset
#'
#' @examples
#' complete_dataset <- generate_complete_dataset(mzrolldb_file_path, metadata)
#' list_complete_dataset <- purrr::map(1:5, ~generate_complete_dataset(mzrolldb_file_path, metadata))
#'
#'@export
generate_complete_dataset <- function(mzroll_list,
                                      metadata,
                                      lod_values) {

  ## merge metadata with mzroll
  ## impute missing peaks
  ## apply median polishing
  mzroll_merged <- claman::merge_samples_tbl(
    mzroll_list = mzroll_list,
    samples_tbl = metadata,
    id_strings = c("SampleName"),
    exact=FALSE) %>%
    impute_missing_peaks(
      lod_values = lod_values,
      quant_var = "log2_abundance") %>%
    claman::normalize_peaks(
      normalization_method = "median polish",
      quant_peak_varname = "log2_abundance",
      norm_peak_varname = "median_log2_abundance")

  ## get group and sample ids of missing peak
  temp_mzroll <- claman::merge_samples_tbl(
    mzroll_list = mzroll_list,
    samples_tbl = metadata,
    id_strings = c("SampleName"),
    exact=FALSE)

  missing_peak_ids <- tidyr::expand_grid(
    groupId = temp_mzroll$features$groupId,
    sampleId = temp_mzroll$samples$sampleId) %>%
    dplyr::anti_join(
      temp_mzroll$measurements,
      by = c("groupId", "sampleId")) %>%
    dplyr::mutate(group_sample_id = paste(groupId, sampleId, sep = "_"))

  ## generate dataframe from mzroll
  output <- as.data.frame(
    romic::triple_to_tidy(mzroll_merged)$data %>%
      ## removing internal standards which have "l" label
      dplyr::filter(!grepl("l", label)) %>%
      dplyr::arrange(ConditionNum) %>%
      dplyr::mutate(condition = factor(condition, levels = unique(condition))) %>%
      dplyr::mutate(SampleName = factor(SampleName, levels = unique(SampleName))) %>%
      dplyr::mutate(group_sample_id = paste(groupId, sampleId, sep = "_")) %>%
      dplyr::mutate(imputed = dplyr::case_when(
        group_sample_id %in% missing_peak_ids$group_sample_id ~ "yes",
        TRUE ~ "no"))
  )
  return(output)
}
