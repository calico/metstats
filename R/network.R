#' Clean a Cytoscape table
#'
#' @param network_name name of the Cytoscape network
#'
#' @return Cytoscape network with the specified table cleaned
#' @export
cytoscape_clean_node_table <- function(network_name = NULL) {

  ## import nodes table to R
  data_nodes <- RCy3::getTableColumns(
    table = "node",
    namespace = "default",
    network = network_name)

  ## check if the required node table columns are present in Cytoscape
  required_columns <- c("KEGG_NODE_TYPE", "KEGG_NODE_LABEL", "KEGG_NODE_LABEL_LIST_FIRST")

  if (!all(required_columns %in% colnames(data_nodes))) {
    stop("Cytoscape node table must have the following columns: KEGG_NODE_TYPE, KEGG_NODE_LABEL, KEGG_NODE_LABEL_LIST_FIRST")
  }

  ## original column names in the Cytoscape node table
  original_columns <- c("SUID", "shared name", "name", "selected", "KEGG_NODE_X",
                        "KEGG_NODE_Y","KEGG_NODE_WIDTH", "KEGG_NODE_HEIGHT",
                        "KEGG_NODE_LABEL", "KEGG_NODE_LABEL_LIST_FIRST",
                        "KEGG_NODE_LABEL_LIST", "id", "KEGG_NODE_LABEL_COLOR",
                        "KEGG_NODE_FILL_COLOR", "KEGG_NODE_REACTIONID",
                        "KEGG_NODE_TYPE", "KEGG_NODE_SHAPE", "KEGG_LINK")

  ## find column names that do not match with the original column names to delete
  columns_to_delete <- colnames(data_nodes)[!colnames(data_nodes) %in% original_columns]

  ## remove previously loaded data from the node table
  if (length(columns_to_delete) > 0) {
    for (i in 1:length(columns_to_delete)) {
      RCy3::deleteTableColumn(columns_to_delete[i], table = "node", namespace = "default", network = network_name)
    }
  }

  data_nodes <- data_nodes %>%
    dplyr::select(tidyselect::matches(original_columns))

  invisible(0)
}

#' Reset Cytoscape network layout to default
#'
#' @description
#' Apply default changes to KEGG pathways in Cytoscape.
#' Since KGML files for KEGG pathways have the same format and properties, same visual modifications can be applied to all KEGG networks in Cytoscape.
#'
#' @param network_name name of the network in Cytoscape
#'
#' @return reset network layout to default
#' @export
cytoscape_mutate_network <- function(network_name = NULL) {

  ## clean node table
  cytoscape_clean_node_table(network_name)

  ## remove node fill colors if they exist
  RCy3::deleteStyleMapping('default', 'NODE_FILL_COLOR')


  ## number of plot slots in Cytoscape is limited to maximum of 9
  plot_slots <- 9

  ## remove previous plots on nodes
  for (i in 1:plot_slots)
  {
    RCy3::removeNodeCustomGraphics(
      slot = i,
      style.name = "default")
  }

  ## select the column for node lables
  RCy3::setNodeLabelMapping(
    'KEGG_NODE_LABEL_LIST_FIRST',
    style.name = "default",
    network = network_name)

  ## change node width
  RCy3::setNodeWidthMapping(
    'KEGG_NODE_TYPE',
    table.column.values = c("ortholog","compound","map"),
    widths = c(5,15,200),
    mapping.type = "d",
    default.width = 1,
    style.name = "default",
    network = network_name)

  ## change node height
  RCy3::setNodeHeightMapping(
    'KEGG_NODE_TYPE',
    table.column.values = c("ortholog","compound","map"),
    heights = c(5,15,30),
    mapping.type = "d",
    default.height = 1,
    style.name = "default",
    network = network_name)

  ## change node font size
  RCy3::setNodeFontSizeMapping(
    'KEGG_NODE_TYPE',
    table.column.values = c("ortholog","compound","map"),
    sizes = c(0,8,15),
    mapping.type = "d",
    default.size = 1,
    style.name = "default",
    network = network_name)

  ## change node shape
  RCy3::setNodeShapeMapping(
    'KEGG_NODE_TYPE',
    table.column.values = c("ortholog","compound","map"),
    c('Rectangle','Circle','Rectangle'),
    default.shape = NULL,
    style.name = "default",
    network = network_name)

  ## modify other visual properties
  RCy3::setNodeColorDefault("#FFFFFF", style.name = "default")
  RCy3::setNodeBorderColorDefault("#000000", style.name = "default")
  RCy3::setNodeBorderWidthDefault(1.5, style.name = "default")
  RCy3::setEdgeColorDefault("#000000", style.name = "default")
  RCy3::setEdgeLineWidthDefault(1, style.name = "default")
  RCy3::setNodeFontFaceDefault("Arial,plain,10", style.name = "default")
  RCy3::setEdgeSourceArrowShapeDefault('ARROW_SHORT', style.name = "default")
  RCy3::setEdgeTargetArrowShapeDefault('ARROW_SHORT', style.name = "default")
  RCy3::setNodeBorderWidthDefault(1, style.name = "default")
  RCy3::setVisualPropertyDefault(list(visualProperty = "NODE_LABEL_POSITION", value = "C,C,c,0,-12"))

  invisible(0)
}

#' Node color visualization
#'
#' @description
#' Visualizing a metric with node color. This functions works well for a comparison of two conditions e.g. log2(treatment/control).
#'
#' @param data_input dataframe of merged Cytoscape table of nodes with experimental data
#' @param vis_metric name of the metric to visualize
#' @param network_name name of the network in Cytoscape
#' @param metric_range range of data for plotting
#'
#' @return Cytoscape network with node color visualizing the metric of interest
#'
#' @export
cytoscape_plot_node_color <- function(data_input,
                                      id_mapping,
                                      vis_metric,
                                      network_name = NULL,
                                      metric_range = 4) {


  ## check if the input data is a dataframe and has a 'feature_name' column
  if (any(class(data_input) == "data.frame")) {
    if (!"feature_name" %in% colnames(data_input)) {
      stop("input data must contain a 'feature_name' column")
    }
  } else {
    stop("input data input must be a data frame")
  }

  all_columns_exist <- all(c("feature_name", "short_name", "id") %in% colnames(id_mapping))


  ## check if the id_mapping is dataframe and has 'feature_name' , 'short_name' and 'id' columns
  if (any(class(id_mapping) == "data.frame")) {
    if (!all(c("feature_name", "short_name", "id") %in% colnames(id_mapping))) {
      stop("id_mapping df must contain 'feature_name', 'short_name', 'id' columns")
    }
  } else {
    stop("id_mapping must be a data frame")
  }

  ## only keep the metric of interest in the data_input
  data_input <- data_input %>%
    dplyr::select("feature_name", vis_metric)

  ## check if the input vis_metric is numeric
  if (class(data_input[[vis_metric]]) != "numeric") {
    stop(paste0("visualization metric must be numeric class\nvis_metric: ",
                vis_metric, "\nclass: ",
                class(data_input[[vis_metric]]))
         )
    }

  ## merge data_input with id_mapping to pull KEGG IDs
  if (!"id" %in% colnames(data_input)) {
    data_input <- dplyr::left_join(data_input %>%
                              dplyr::mutate(feature_name = tolower(feature_name)),
                            id_mapping %>% dplyr::select(feature_name, short_name, id) %>%
                              dplyr::mutate(feature_name = tolower(feature_name)), by = "feature_name")
    }

  ## reset network layout to default
  cytoscape_mutate_network(network_name)

  ## import nodes table from Cytoscape to R
  data_nodes <- RCy3::getTableColumns(
    table = "node",
    namespace = "default",
    network = network_name)

  ## add 'exists' column to distinguish measured vs. missing metabolites
  data_input$exists <- data_input$id

  ## merge experimental data with nodes table
  data_nodes_mut <- dplyr::left_join(data_nodes, data_input, by = c("KEGG_NODE_LABEL" = "id"))
  data_nodes_mut$exists <- ifelse(!is.na(data_nodes_mut$exists), "yes", "no")

  data_nodes_mut <- data_nodes_mut %>%
    dplyr::mutate(short_name = dplyr::case_when(
      is.na(short_name) ~ KEGG_NODE_LABEL_LIST_FIRST,
      short_name == "" ~ KEGG_NODE_LABEL_LIST_FIRST,
      TRUE==TRUE ~ short_name))

  col_identifier <- match("feature_name", names(data_nodes_mut)) + 1
  data_nodes_mut[,(col_identifier:ncol(data_nodes_mut))][is.na(data_nodes_mut[,(col_identifier:ncol(data_nodes_mut))])] = 0


  ## change KEGG_NODE_TYPE for missing compounds
  data_nodes_mut <- data_nodes_mut %>%
    dplyr::mutate(node_type = dplyr::case_when(
      grepl("no", exists) & grepl("compound", KEGG_NODE_TYPE) ~ "missing",
      TRUE==TRUE ~ KEGG_NODE_TYPE
    ))

  ## add * to the name of imputed compounds
  data_nodes_mut$node_name <- data_nodes_mut$short_name
  if (!is.na(match("imputed",names(data_nodes_mut)))) {
    data_nodes_mut <- data_nodes_mut %>%
      dplyr::mutate(node_name = ifelse(imputed == "yes", paste(node_name, "*", sep = " "), node_name))
    }

  ## add legend for the plot
  data_nodes_mut$node_name[data_nodes_mut$KEGG_NODE_LABEL == "Legend"] <- vis_metric


  ## load experimental data to Cytoscape
  RCy3::loadTableData(
    data_nodes_mut,
    data.key.column = "KEGG_NODE_LABEL",
    table = "node",
    table.key.column = "KEGG_NODE_LABEL",
    namespace = "default",
    network = network_name)

  ## select the column for node lables
  RCy3::setNodeLabelMapping(
    'node_name',
    style.name = "default",
    network = network_name)

  ## change node font size
  RCy3::setNodeFontSizeMapping(
    'node_type',
    table.column.values = c("ortholog","compound","missing","map"),
    sizes = c(0,7,7,15),
    mapping.type = "d",
    default.size = 1,
    style.name = "default",
    network = network_name)

  ## function to visualize a metric with node color
  RCy3::setNodeColorMapping(
    vis_metric,
    table.column.values = c(-metric_range, 0, metric_range),
    colors = c('#0000FF','#FFFFFF','#FF3333'), #blue-white-red color scale for negative-zero-positive values
    mapping.type = "c",
    default.color = NULL,
    style.name = "default",
    network = network_name)


  ## change node shape
  RCy3::setNodeShapeMapping(
    'node_type',
    table.column.values = c("ortholog","compound","missing","map"),
    c('Rectangle','Circle','Diamond','Rectangle'),
    default.shape = NULL,
    style.name = "default",
    network = network_name)

  ## change node width
  RCy3::setNodeWidthMapping(
    'node_type',
    table.column.values = c("ortholog","compound","missing","map"),
    widths = c(5,15,15,200),
    mapping.type = "d",
    default.width = 1,
    style.name = "default",
    network = network_name)

  ## change node height
  RCy3::setNodeHeightMapping(
    'node_type',
    table.column.values = c("ortholog","compound","missing","map"),
    heights = c(5,15,15,30),
    mapping.type = "d",
    default.height = 1,
    style.name = "default",
    network = network_name)
}
