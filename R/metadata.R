#' Import metadata
#'
#' @returns
#'   a metadata table
#'
#' @export
import_metadata <- function(sheet_id) {

  metadata_tbl <- googlesheets4::read_sheet(
    sheet_id,
    sheet = "metadata"
    ) %>%
    dplyr::mutate_if(is.list, as.character)

  ## validate metadata
  check_metadata(metadata_tbl)

  return(metadata_tbl)
  }

#' Check metadata
#'
#' @returns
#'   nothing, only checks
#' @export
check_metadata <- function(metadata_tbl) {

  ## require a few standard fields in metadata
  required_sample_vars <- c(
    "SampleName",
    "Replicate",
    "ConditionNum",
    "ReferenceConditionNum"
    )

  missing_required_vars <- setdiff(
    required_sample_vars,
    colnames(metadata_tbl)
    )

  if (length(missing_required_vars) != 0) {
    stop(
      "\"Sample List\" is missing required variables: ",
      paste(missing_required_vars, collapse = ", ")
      )
    }

  ## check for SampleName uniqueness
  duplicated_SampleName <- metadata_tbl %>%
    dplyr::count(SampleName) %>%
    dplyr::filter(n > 1)

  if (nrow(duplicated_SampleName) != 0) {
    stop(
      nrow(duplicated_SampleName),
      "SampleName were not unique, duplicated names: ",
      paste(duplicated_SampleName$SampleName, collapse = ", ")
    )
    }

  ## check that all reference conditions exist
  if (
    !all(metadata_tbl$ReferenceConditionNum %in% metadata_tbl$ConditionNum)
    ) {
    stop("some ReferenceConditionNum not found as ConditionNum")
    }

  return(invisible(0))
  }
