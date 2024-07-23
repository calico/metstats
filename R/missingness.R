#' Missing value imputation
#'
#'@description impute missing peaks with provided feature-level imputation values
#'
#'@param mzroll_list: data in triple omic structure
#'@param lod_values: a tibble that maps groupId to feature-level imputation values
#' if a tibble is not provided, half min value per feature will be used as imputation values
#'@param quant_var: column to use for peak values, the function performs imputation in log2 space, peak quant values must be in log2 space
#'@param imputation_sd: standard deviation of Gaussian distribution to use for missing peak imputation
#'
#'@return triple omic data with imputed missing peaks
#'
#'@export
impute_missing_peaks <- function(mzroll_list,
                                 lod_values = NULL,
                                 quant_var = "log2_abundance",
                                 imputation_sd = 0.15) {

  valid_quant_var <- dplyr::setdiff(
    mzroll_list$design$measurements$variable,
    c(mzroll_list$design$feature_pk, mzroll_list$design$sample_pk)
  )

  checkmate::assertChoice(quant_var, valid_quant_var)
  checkmate::assertNumeric(mzroll_list$measurements[[quant_var]])

  ##  get half min value per feature to use for imputation if feature-specific imputation values are not provided
  if(is.null(lod_values)) {
    lod_values <- mzroll_list$measurements %>%
      dplyr::group_by(groupId) %>%
      dplyr::filter(!!rlang::sym(quant_var) == min(!!rlang::sym(quant_var))/2) %>%
      dplyr::ungroup() %>%
      dplyr::select(groupId, !!rlang::sym(quant_var))
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
    dplyr::mutate(!!rlang::sym(quant_var) := stats::rnorm(1, mean = !!rlang::sym(quant_var)+1, sd = imputation_sd))

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
#'@param feature_id groupId from mzrolldb
#'@param cond Condition in metadata
#'@param metadata metadata of experiment
#'@param df dataframe in long format
#'
#'@return dataframe of terms and groupIds for comparisons with imputed values
#'
#'@export
imputed_comparisons <- function(feature_id,
                                cond,
                                metadata,
                                df) {

  output <- data.frame(
    term = character(1),
    groupId = factor(1),
    imputed = "no")

  if (nrow(metadata %>% filter(Condition == cond)) < 2) {return(output)} else {
    ref_cond <- unique(metadata$RefCondition[metadata$Condition == cond])
    cond <- unique(metadata$Condition[metadata$Condition == cond])

    if (cond == ref_cond) {return(output)} else {
      df <- df %>%
        dplyr::filter(Condition %in% c(cond, ref_cond)) %>%
        dplyr::mutate(Condition = factor(Condition, levels = c(ref_cond, cond))) %>%
        dplyr::filter(groupId == feature_id) %>%
        dplyr::select(imputed)
      if ("yes" %in% df$imputed) {
        output$term = paste("Condition", cond, sep = "")
        output$groupId = feature_id
        output$imputed = "yes"
      }
    }
    return(output)
    }
  }

#' Generate complete dataset upon imputing missing values
#'
#'@param mzroll_list: data in triple omic structure that's processed with claman::process_mzroll, samples table of mzroll_list must contain a numeric InjOrder column
#'@param metadata: metadata of experiment
#'@param lod_values: a tibble that maps groupId to log2 feature-level imputation values
#'
#'@return dataframe of complete dataset
#'
#'@export
generate_complete_dataset <- function(mzroll_list,
                                      metadata,
                                      lod_values) {

  ## check that InjOrder column exists in mzroll_list and is numeric
  stopifnot("InjOrder" %in% colnames(mzroll_list$samples))
  checkmate::assertNumeric(mzroll_list$samples[["InjOrder"]])

  ## check that SampleName column exists in metadata
  stopifnot("SampleName" %in% colnames(metadata))

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
      normalization_method = "lm",
      quant_peak_varname = "log2_abundance",
      norm_peak_varname = "log2_abundance_lm",
      time_col_varname = "InjOrder") %>%
    claman::normalize_peaks(
      normalization_method = "median polish",
      quant_peak_varname = "log2_abundance_lm",
      norm_peak_varname = "log2_abundance_lm_median") %>%
    claman::median_polish_predict_dilutions() %>%
    claman::normalize_peaks(
      normalization_method = "center",
      quant_peak_varname = "log2_abundance_lm_median",
      norm_peak_varname = "log2_normalized")

  mzroll_merged$measurements$log2_normalized <- as.numeric(mzroll_merged$measurements$log2_normalized)

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
      dplyr::arrange(Condition) %>%
      dplyr::mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
      dplyr::mutate(SampleName = factor(SampleName, levels = unique(SampleName))) %>%
      dplyr::mutate(group_sample_id = paste(groupId, sampleId, sep = "_")) %>%
      dplyr::mutate(imputed = dplyr::case_when(
        group_sample_id %in% missing_peak_ids$group_sample_id ~ "yes",
        TRUE ~ "no"))
  )
  return(output)
}
