#'Pathway Enrichment Analysis
#'
#'@param data_input:
#'      dataframe of experimental data in wide format, must contain "feature_name" and "ranking_metric" columns
#'
#'@param id_mapping:
#'      dataframe with matching compound names and ids (such as kegg ids), ids should match with ids in feature_sets
#'      must contain "feature_name" and "id" columns
#'
#'@param feature_sets
#'      path to the feature sets file that should have a standard GSEA .gmt file format, see the link below:
#'      https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
#'      for MSEA, this file is embedded within the package:
#'      "/metstats/inst/kegg_msea.txt"
#'
#'@param ranking_metric
#'      column name from input data to be used as the metric for enrichment analysis
#'
#'@param metric_type
#'      type of the metric being used for enrichment analysis "abs" or "non-abs"
#'
#'@param pval
#'      pvalue threshold [DEFAULT = 1]
#'
#'@description
#'       Pathway Enrichment Analysis
#'
#'@export
analysis_pathway_enrichment <- function(data_input,
                                        id_mapping,
                                        feature_sets = NULL,
                                        ranking_metric,
                                        metric_type = "non-abs",
                                        pval = 1) {

  ## TODO: figure out where to put the .txt file within the repo
  if (is.null(feature_sets)) {
    #feature_sets = system.file("kegg_msea_set_revised.txt", package="metstats")
    feature_sets = "kegg_msea.txt"
  }

  ## check if the input data is a df and has a 'feature_name' column
  if (any(class(data_input) == "data.frame")) {
    if (!"feature_name" %in% colnames(data_input)) {
      stop("input must contain a 'feature_name' column")
    }
  } else {
    stop("input must be a data frame")
  }

  ## map experimental data with compound ids
  if (!"id" %in% colnames(data_input)) {
    data_input <- dplyr::left_join(data_input %>%
                              dplyr::mutate(feature_name = tolower(feature_name)),
                            id_mapping %>%
                              dplyr::select(feature_name,id) %>%
                              dplyr::mutate(feature_name = tolower(feature_name)),
                            by = "feature_name")
    }

  ## remove compounds with missing id

  data_input <- data_input %>%
    dplyr::filter(!(is.na(data_input$id) | data_input$id == ""))

  ranked_list <- data_input %>%
    dplyr::mutate(ranking_metric := !! dplyr::sym(ranking_metric)) %>%
    dplyr::select(id, ranking_metric)

  names(ranked_list)[names(ranked_list) == ranking_metric] <- "ranking_metric"
  ranked_mod = ranked_list$ranking_metric
  names(ranked_mod) = ranked_list$id
  ranked_mod = sort(ranked_mod, decreasing = TRUE)
  ranked_mod = ranked_mod[!duplicated(names(ranked_mod))]

  if ( any( duplicated(names(ranked_mod)) )  ) {
    warning("Duplicates in the list of compounds")
    ranked_mod = ranked_mod[!duplicated(names(ranked_mod))]
  }

  if  ( !all( order(ranked_mod, decreasing = TRUE) == 1:length(ranked_mod)) ) {
    ranked_mod = sort(ranked_mod, decreasing = TRUE)
  }

  enrich_sets = fgsea::gmtPathways(feature_sets)

  fgRes <- fgsea::fgsea(pathways = enrich_sets,
                        stats = ranked_mod) %>%
    as.data.frame() %>%
    dplyr::filter(padj < !!pval)

  fgRes = fgRes %>% dplyr::arrange(dplyr::desc(NES))

  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")

  ## select top 10 up and top 10 down sets
  ## or select top 20 sets if the ranking metric being used is absolute
  if (metric_type == "abs") {
    filtRes = fgRes[fgRes$NES > 0,]
    filtRes = utils::head(filtRes, n = 20)
  } else {
    filtRes = rbind(utils::head(fgRes, n = 10),
                    utils::tail(fgRes, n = 10 ))
  }


  filtRes <- filtRes %>% dplyr::arrange(NES)
  filtRes$pathway <- factor(filtRes$pathway, levels=unique(filtRes$pathway))

  table_unnested <- fgRes %>%
    tidyr::unnest(cols = leadingEdge)

  compound_name_id <- id_mapping %>%
    dplyr::select(feature_name, id) %>%
    unique()

  table_updated <- table_unnested %>%
    dplyr::left_join(compound_name_id, by =c("leadingEdge"="id")) %>%
    dplyr::mutate(name_or_id = ifelse(is.na(feature_name), leadingEdge, feature_name)) %>%
    dplyr::select(-leadingEdge, -feature_name) %>%
    dplyr::rename(leadingEdge = name_or_id)

  leadingEdge <- stats::aggregate(leadingEdge ~ pathway, table_updated, paste)

  fgRes_modified <- fgRes %>%
    dplyr::select(-leadingEdge)

  if (metric_type == "abs") {
    fgRes_modified = fgRes_modified[fgRes_modified$NES > 0,]
  }

  enrichment_table <- dplyr::left_join(fgRes_modified, leadingEdge, by = "pathway")

  ## add size of pathway to the enrichment table
  pathway_size <- as.data.frame(names(enrich_sets)) %>%
    dplyr::rename("pathway" = "names(enrich_sets)") %>%
    dplyr::mutate(setsize = NA)

  for (i in 1:nrow(pathway_size)) {
    pathway_size[i,2] <- length(enrich_sets[[as.character(pathway_size[i,1])]][enrich_sets[[as.character(pathway_size[i,1])]]!=""])
  }

  enrichment_table <- dplyr::left_join(enrichment_table, pathway_size, by = "pathway") %>%
    dplyr::relocate(setsize, .after = size)

  ## plot a dot plot of all enriched sets
  if (metric_type == "abs") {
    enrichment_plot <- ggplot2::ggplot(filtRes, aes(x = NES, y = pathway, size = -log10(pval))) +
      ggplot2::geom_point(aes(colour = NES, size = -log10(pval))) +
      ggplot2::scale_color_gradient2(low="white",
                            high="red", space ="Lab", limits = c(0, plyr::round_any(max(abs(filtRes$NES)),0.5,f = ceiling))) +
      labs(y="Pathway", x="Normalized Enrichment Score", title=ranking_metric) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1.5), panel.background = ggplot2::element_blank(),
            axis.line = element_line(colour = "black", size = 0.65),
            axis.text.x = ggplot2::element_text(size=15, angle = 0, colour = "black", face = "bold"),
            axis.text.y = ggplot2::element_text(size=15,  colour = "black", face = "bold"),
            axis.ticks.x = element_line(colour="black"),
            axis.ticks.y = element_line(colour="black"),
            axis.title.y = ggplot2::element_text(size=20, face = "bold", colour="black"),
            axis.title.x = ggplot2::element_text(size=20, face = "bold", colour="black"),
            legend.text=ggplot2::element_text(size=12, face = "bold"),
            plot.title = ggplot2::element_text(size=20, face = "bold"),
            strip.text.x = ggplot2::element_text(size = 15, face = "bold"))
    } else {
    enrichment_plot <- ggplot2::ggplot(filtRes, aes(x = NES, y = pathway, size = -log10(pval))) +
      ggplot2::geom_point(aes(colour = NES, size = -log10(pval))) +
      ggplot2::scale_color_gradient2(low="blue", mid="white",
                            high="red", space ="Lab", limits = c(-plyr::round_any(max(abs(filtRes$NES)),0.5,f = ceiling), plyr::round_any(max(abs(filtRes$NES)),0.5,f = ceiling))) +
      ggplot2::labs(y="Pathway", x="Normalized Enrichment Score", title=ranking_metric) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1.5), panel.background = ggplot2::element_blank(),
            axis.line = element_line(colour = "black", size = 0.65),
            axis.text.x = ggplot2::element_text(size=15, angle = 0, colour = "black", face = "bold"),
            axis.text.y = ggplot2::element_text(size=15,  colour = "black", face = "bold"),
            axis.ticks.x = element_line(colour="black"),
            axis.ticks.y = element_line(colour="black"),
            axis.title.y = ggplot2::element_text(size=20, face = "bold", colour="black"),
            axis.title.x = ggplot2::element_text(size=20, face = "bold", colour="black"),
            legend.text=ggplot2::element_text(size=12, face = "bold"),
            plot.title = ggplot2::element_text(size=20, face = "bold"),
            strip.text.x = ggplot2::element_text(size = 15, face = "bold"))
    }


  # generate mountain plots of enriched sets
  mountain_plot <- list()
  for (i in 1:length(fgRes$pathway)) {
    mountain_plot[[i]] <- fgsea::plotEnrichment(enrich_sets[[fgRes$pathway[i]]],
                                                ranked_mod) + ggplot2::ggtitle(fgRes$pathway[i], ranking_metric)
  }

  names(mountain_plot) <- fgRes$pathway

  output = list("enrichment_table" = enrichment_table, "enrichment_plot" = enrichment_plot, "mountain_plot" = mountain_plot)
  return(output)
}
