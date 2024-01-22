#' Table of peakgroupmatch from file
#'
#'@param mzrolldb_file_path: file path to mzrolldb file
#'
#'@return compounds from mzrolldb file
#'
#'@export
PDB_compounds <- function(mzrolldb_file_path) {
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = mzrolldb_file_path)

  samples <- dplyr::tbl(con, "compounds") %>%
    collect() %>%
    dplyr::arrange(compoundId)

  DBI::dbDisconnect(conn = con)

  return(samples)
}
