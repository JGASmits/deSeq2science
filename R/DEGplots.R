#' GOenr_simplified
#'
#' This function runs GO enrichment using clusterprofilea nd simplify on the outcome
#' @param gene_list list of genes (DEGS)
#' @param output_file name of the output (pdf) file
#' @param res DESeq2 result of differential testing, used for selecting most significant genes from gene_list for plotting
#' @param dds deseq2 object, used for grabbing count data of genes of GO-term
#' @param fill_annotation metadata column deseq2 object used for count plotting fill
#' @param x_axis_annotation metadata column deseq2 object used for count seperation x-axis
#' @param topn_genes n genes plotted per GO-term
#' @param topn_GOterms n GO-terms plotted
#' @param gene_background_list optional gene background
#' @param simplify_value value wich is used to simplify go terms, smaller number mean less GO-concationation.
#' @export
GOenr_simplified <- function(gene_list,
                              output_file,
                              res,
                              dds,
                              fill_annotation = 'disease',
                              x_axis_annotation = 'timepoint',
                              topn_genes = 9,
                              topn_GOterms = 20,
                              gene_background_list = 'NA',
                              qval_cutoff = 0.05,
                              gene_ID = 'SYMBOL',
                              GO_ontology="BP",
                              simplify_value = 0.7 ){
  #run clusterprofiler GO-term enrichment on the genelist
  #simplify's the GO term results
  #return a list with 1. a GO_term object, 2 a df of the GO_term object
  #write GO_term output table to a csv file in output_file
  #library(clusterProfiler)
  #library(org.Hs.eg.db)
  #library(tidyverse)
  if (gene_background_list == "NA"){
    GO_obj <- clusterProfiler::enrichGO(gene = gene_list,
                                        OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                        keyType       = gene_ID,
                                        ont           = GO_ontology,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff  = qval_cutoff)
  }else{
    GO_obj <- clusterProfiler::enrichGO(gene = gene_list,
                                        universe = gene_background_list,
                                        OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                        keyType       = gene_ID,
                                        ont           = GO_ontology,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff  = qval_cutoff)}

  if(nrow(as.data.frame(GO_obj))==0){
    return('no sig go terms')}else{
      GO_obj_s <- clusterProfiler::simplify(GO_obj, cutoff = simplify_value)
      GO_df <- as.data.frame(GO_obj_s)
      #stringr::str_split(GO_df$GeneRatio, '/')
      GO_df <- GO_df %>% tidyr::separate(GeneRatio, c("genes_found", "GO_size"), "/")
      GO_df$genes_found <- as.numeric(GO_df$genes_found)
      GO_df <- GO_df[order(GO_df$genes_found, decreasing = T),]
      write.table(as.data.frame(GO_df), file= paste0(output_file,'.csv'), sep = ',', row.names = T)

      pdf(paste0(paste0(output_file,'.pdf')) ,width=12,height=12,paper='special')
      print(clusterProfiler::dotplot(GO_obj_s,showCategory = topn_GOterms))#print overview dotplot

      if(dim(GO_df)[1] > topn_GOterms){max_examples <- topn_GOterms}else{max_examples <- dim(GO_df)[1]}
      for (i in 1:max_examples){
        count_plot <- NULL
        genes_enriched <- stringr::str_split(GO_df[i,]$geneID,'/')[[1]]
        genes_enriched_res <- res[genes_enriched,]
        top <- rownames(head(genes_enriched_res[order(genes_enriched_res$padj),],topn_genes))  #take top X  differential genes associated to the GO_term
        tcounts <- t(log2((counts(dds[top, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
          merge(colData(dds), ., by="row.names") %>%
          tidyr::gather(gene, expression, (ncol(.)-length(top)+1):ncol(.))

        count_plot <- ggplot2::ggplot(tcounts, ggplot2::aes_string(x_axis_annotation, 'expression', fill = fill_annotation)) +
          ggplot2::geom_boxplot() +
          ggplot2::facet_wrap(~gene, scales="free_y") +
          ggplot2::labs(x=GO_df[i,]$Description,
                        y="Expression (log normalized counts)",
                        title=paste0("Top GO-term genes ",GO_df[i,]$Description))

        print(count_plot)}
      dev.off()
    }
}
