#' Plot BMD, BMDL, and BMDU for a Set of Genes
#'
#' This function generates a plot showing BMD, BMDL, and BMDU values for a given set of genes.
#'
#' @param BMDFilMat A data frame containing BMD values for multiple genes.
#' @param gi A vector of gene names for which to plot BMD values.
#'
#' @return A plot showing BMD, BMDL, and BMDU values for the selected set of genes.
#' @export
plot_bmd_bmdl_bmdu_set_of_genes = function(BMDFilMat,
                                           gi,
                                           enrich_ppi_info = TRUE,
                                           gene_id_type = "ENSEMBL",
                                           organism = "human"){

  idx = which(tolower(BMDFilMat[,"Feature"]) %in% tolower(gi))
  BMD = as.numeric(as.vector(BMDFilMat[idx,"BMD"]))
  BMDL = as.numeric(as.vector(BMDFilMat[idx,"BMDL"]))
  BMDU = as.numeric(as.vector(BMDFilMat[idx,"BMDU"]))

  names(BMD) = BMDFilMat[idx,"Feature"]
  names(BMDL) = names(BMDU) = names(BMD)


  BMD = data.frame(gene = names(BMD), bmd = BMD, bmdl = BMDL, bmdu = BMDU)
  to_rem = which(is.na(BMD$bmd) | is.na(BMD$bmdl) | is.na(BMD$bmdu))
  if (length(to_rem) > 0) BMD = BMD[-to_rem,]
  BMD = BMD[order(BMD$bmd),]

  BMD = BMD %>%  dplyr::group_by(gene) %>% dplyr::summarise( dplyr::across(everything(), list(mean)))
  colnames(BMD) = c("gene","bmd","bmdl","bmdu")
  BMD$gene = factor(x = BMD$gene, levels = BMD$gene)
  BMD$gene <- forcats::fct_reorder(BMD$gene, BMD$bmd, .desc = FALSE)
  BMD= as.data.frame(BMD)
  rownames(BMD)= BMD$gene
  BMD=BMD[BMD$gene,]

  if(enrich_ppi_info){
    library(igraph)
    ppi_g = load_ppi_data(organism,gene_id_type,graph_type = "ppi")
    tf_g = load_ppi_data(organism,gene_id_type,graph_type = "tf")
    mirna_g = load_ppi_data(organism,gene_id_type,graph_type = "mirna")

    degree_nodes = igraph::degree(ppi_g)

    tf_list    <- unique(ends(tf_g, E(tf_g), names = T)[, 1])
    mirna_list <- unique(ends(mirna_g, E(mirna_g), names = FALSE)[, 1])

    BMD = cbind(BMD, "degree" = degree_nodes[rownames(BMD)],
                "is_TF" = rownames(BMD) %in% tf_list,
                "is_mirna" = rownames(BMD) %in% mirna_list)

    # p = ggplot2::ggplot(data = as.data.frame(BMD), ggplot2::aes(x = log2(degree), y = bmd, group = 1, label1 = bmdl, label2 = bmdu)) +
    # ggplot2::geom_line() +
    # ggplot2::geom_ribbon(aes(ymin=bmdl, ymax = bmdu), linetype = 2, alpha = 0.1) +

    p = ggplot2::ggplot(data = as.data.frame(BMD),
                        ggplot2::aes(x = log2(degree), y = bmd,
                                     color = is_TF)) +
      ggplot2::geom_point() +
      ggplot2::geom_text(ggplot2::aes(label = gene),
                         hjust = -0.1, vjust = 0, size = 3) +
      ggplot2::labs(y = "BMD", x = "log2 degree PPI") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, vjust = 1, hjust = 1))

  }else{
    p = ggplot2::ggplot(data = as.data.frame(BMD), ggplot2::aes(x = gene, y = bmd, group = 1, label1 = bmdl, label2 = bmdu)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_ribbon(aes(ymin=bmdl, ymax = bmdu), linetype = 2, alpha = 0.1) +
      ggplot2::labs(y = "BMDL - BMD - BMDU", x = "Gene") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, vjust = 1, hjust = 1))
  }


  return(p)
}
