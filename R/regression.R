#' Pairwise linear regression for one feature and one reference condition
#' This will perform regression for all conditions that have the same reference condition simultaneously
#'
#'@param feature_id: groupId from mzrolldb
#'@param cond_num: ConditionNum in metadata
#'@param metadata: metadata of experiment
#'@param df: dataframe in long format
#'
#'@return one-row data from linear regression
#'
#' @examples
#' complete_dataset <- generate_complete_dataset(mzrolldb_file_path, metadata)
#' list_complete_dataset <- purrr::map(1:5, ~generate_complete_dataset(mzrolldb_file_path, metadata))
#'
#'@export
lm_feature <- function(feature_id, ref_cond_num, metadata, df) {

  cond_nums <- metadata %>%
    dplyr::filter(ReferenceConditionNum == ref_cond_num & ConditionNum != ref_cond_num) %>%
    dplyr::select(ConditionNum) %>%
    dplyr::group_by(ConditionNum) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    unique() %>%
    unlist() %>%
    as.numeric()

  ## skip regression analysis if there is less than two sample left
  if (nrow(metadata %>% dplyr::filter(ConditionNum %in% cond_nums)) < 2) {return(NULL)} else {

    ref_cond <- metadata %>%
      dplyr::filter(ConditionNum == ref_cond_num) %>%
      dplyr::select(condition) %>%
      unique() %>%
      unlist() %>%
      as.character()

    conds <- metadata %>%
      dplyr::filter(ConditionNum %in% cond_nums) %>%
      dplyr::select(condition) %>%
      unique() %>%
      unlist() %>%
      as.character()

    ## subset data
    data_temp <- df %>%
      dplyr::filter(ConditionNum %in% c(cond_nums, ref_cond_num)) %>%
      dplyr::mutate(condition = factor(condition, levels = c(
        metadata %>%
          dplyr::filter(ConditionNum == ref_cond_num) %>%
          dplyr::select(condition) %>%
          unique() %>%
          unlist() %>%
          as.character(),
        metadata %>%
          dplyr::filter(ConditionNum %in% cond_nums) %>%
          dplyr::select(condition) %>%
          unique() %>%
          unlist() %>%
          as.character()
      )
      )
      ) %>%
      dplyr::filter(groupId == feature_id)

    ## linear regression
    output <- with(data_temp , stats::lm(log2_abundance ~ condition))
    return(output)
  }
}

#' Linear regression for multiple conditions
#'
#'@param feature_id: groupId(s) from mzrolldb
#'@param ref_condition_nums: RefConditionNum(s) in metadata
#'@param metadata: df of metadata
#'@param df: data in long format
#'
#'@return one-row data from linear regression
#'
#' @examples
#' lm_all <- purrr::map(feature_ids, ~lm_all(.x, ref_ids, metadata, complete_data))
#'
#'@export
lm_multi <- function(feature_id, ref_condition_nums, metadata, df) {
  output <- data.frame()
  for (i in ref_condition_nums) {
    lm_list <- lm_feature(feature_id, i, metadata, df)
    if(!is.null(lm_list)) {
      output <- rbind(
        output,
        as.data.frame(coef(summary(lm_list))) %>%
          dplyr::mutate(term = row.names(.)) %>%
          dplyr::relocate(term) %>%
          dplyr::filter(term !="(Intercept)") %>%
          dplyr::mutate(groupId = feature_id)
      )
    }
  }
  colnames(output) <- c("term", "estimate", "std.error", "statistics", "p.value", "groupId")
  return(output)
}

#' Pooled results for linear regression of one feature and all pairwise comparisons upon multiple imputation of missing peaks
#'
#'@param feature_id: groupId(s) from mzrolldb
#'@param ref_condition_nums: RefConditionNum(s) in metadata
#'@param metadata: metadata of experiment
#'@param df_list: list of complete dataframes upon imputing missing peaks
#'
#'@return dataframe of linear regression summary for one feature and all pairwise comparisons
#'
#' @examples
#' lm_pooled <- lm_pool(2, c(1,2,3), metadata, df_list)
#'
#'@export
lm_pool <- function(feature_id, ref_condition_nums, metadata, df_list) {
  lm_temp <- data.frame(
    term = character(),
    estimate = numeric(),
    std.error = numeric(),
    statistics = numeric(),
    df = numeric(),
    p.value = numeric())
  for (i in ref_condition_nums) {
    imputation_list <- purrr::map(df_list, ~lm_feature(feature_id, i, metadata, .x))

    if(!is.null(imputation_list[[1]])) {
      lm_temp <- rbind(
        lm_temp,
        summary(mice::pool(imputation_list)) %>%
          dplyr::filter(term !="(Intercept)") %>%
          dplyr::mutate(groupId = feature_id)
      )
    }
  }
  return(lm_temp)
}

#' Calculate qvalues
#'
#'@param term_data: dataframe of linear regression conataining p.value
#'
#'@return dataframe of linear regression with qvalues added
#'
#' @examples
#' lm_fdr <- fdr(lm_data)
#'
#'@export
fdr <- function(term_data) {
  p_values <- term_data$p.value
  q_values <- try(qvalue::qvalue(p_values)$qvalues, silent = TRUE)

  if ("try-error" %in% class(q_values)) {
    # if qvalue fails this is probably because there are no p-values greater
    # than 0.95 (the highest lambda value)
    # if so add a single p-value of 1 to try to combat the problem
    q_values <- qvalue::qvalue(c(p_values, 1))$qvalues
    q_values <- q_values[-length(q_values)]
  }

  term_data %>%
    dplyr::mutate(qvalue = q_values)
}

#' Volcano plot
#'
#' @param lm_data: dataframe of regression analysis
#' @param FDR_cutoff: FDR cutoff to label for significance
#' @param feature_labels: list of compound names from compoundName column to label on the plot
#'
#' @returns a grob
#'
#' @examples
#' volcano_plot(lm_data, 0.05, c("Serine", "Glycine"))
#'
#' @export
volcano_plot <- function(
  lm_data,
  FDR_cutoff = 0.05,
  feature_labels = NULL
) {

  checkmate::assertDataFrame(lm_data)
  stopifnot("term" %in% colnames(lm_data))

  effect_var <- dplyr::case_when(
    "estimate" %in% colnames(lm_data) ~ "estimate",
    "sumsq" %in% colnames(lm_data) ~ "sumsq",
    TRUE ~ NA_character_
  )

  if (is.na(effect_var)) {
    stop("volcano plot cannot be generated due to unknown test")
  }

  lm_data %>%
    dplyr::filter(!is.na(p.value)) %>%
    dplyr::mutate(
      p.value.trans = -log10(p.value),
      is_discovery = qvalue < FDR_cutoff,
      plot_color = dplyr::case_when(
        imputed == "yes" ~ "imputed",
        is_discovery == TRUE ~ "significant",
        TRUE ~ "ns"
      )
    ) %>%
    ggplot(aes_string(x = effect_var)) +
    {if ("compoundName" %in% colnames(lm_data)) {
      geom_point(aes(y = p.value.trans, color = plot_color, name = compoundName))
    }
      else {
        geom_point(aes(y = p.value.trans, color = plot_color))
      }
    } +
    {if ("compoundName" %in% colnames(lm_data)) {
      geom_text(aes(label = ifelse(compoundName %in% feature_labels, compoundName, ""), y = p.value.trans, vjust = -0.75))
    }
    } +
    scale_x_continuous("Effect size") +
    scale_y_continuous(expression(-log[10] ~ "pvalue")) +
    scale_color_manual(values = c("ns" = "gray50", "significant" = "red", "imputed" = "blue")) +
    theme_bw()
}
