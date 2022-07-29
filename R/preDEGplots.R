#' plot sample distances
#'
#' This function plots the sample distance of samples within a deseq2 object
#' @param normalized_deseq_int rlog or vst normalized deseq2 object
#' @param output_file name of the output (pdf) file
#' @param annotation1 potential name of metadata1 in the dds object plotted as annotation
#' @param annotation2 potential name of metadata2 in the dds object plotted as annotation
#' @param costum_sample_order list of all sample names in the desired order
#' @export
plot_dist_dds <- function(normalized_deseq_int, output_file, annotation1 = NA, annotation2 = NA, flip_tree_coordinates= NA, costum_sample_order = NA) {
 #output_file = location in which the output pdf is generated
  #annotation1 = optional annotation name, should be present as metadata to the normalized_deseq_int object
  #annotation2 = optional second annotation name, should be present as metadata to the normalized_deseq_int object
  sampleDists <- dist(t(SummarizedExperiment::assay(normalized_deseq_int)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9,"Blues")) )(1000)
  
  dend <- stats::as.dendrogram(stats::hclust(as.dist(sampleDistMatrix)))
  #dend = dend %>% dendextend::set("labels_to_char") %>% labels
  
  if (typeof(costum_sample_order) == "character"){
    print('reordering dendogram based on costum sample order list')
    dend <- dend %>% dendextend::rotate(costum_sample_order)} #rotate to match labels new order}
  
  if (is.na(annotation1)) {
    myplot <- 
      ComplexHeatmap::Heatmap(sampleDistMatrix,
                              row_dend_side = "right",
                              show_column_dend = FALSE,
                              cluster_rows = dend,
                              cluster_columns = dend,
                              col = colors)}else {
                                if (is.na(annotation2)){
                                  ha = ComplexHeatmap::rowAnnotation(annotation1 = normalized_deseq_int[[annotation1]])}
                                else {ha = ComplexHeatmap::rowAnnotation(annotation1 = normalized_deseq_int[[annotation1]],
                                                                         annotation2 = normalized_deseq_int[[annotation2]])}
                                myplot <- 
                                  ComplexHeatmap::Heatmap(sampleDistMatrix,
                                                          row_dend_side = "right",
                                                          show_column_dend = FALSE,
                                                          col = colors,
                                                          cluster_rows = dend,
                                                          cluster_columns = dend,
                                                          right_annotation = ha )}
  
  pdf(output_file,width=6,height=6,paper='special') 
  print(myplot)
  dev.off()
}
#' plot PCA of samples
#'
#' This function plots the sample distance of samples within a deseq2 object
#' @param object deseq2 object
#' @param intgroup potential metadata columns, e.g. c('genotype', 'timepoint')
#' @param output_dir directory where pdf files will be generated
#' @param run_go_enrichment True or False wether to run GO_enrichment using clusterprofiler
#' @param Genenames_keyType genename type, currently only 'SYMBOL'is supported
#' @param run_go_enrichment True or False wether to run GO_enrichment
#' @param go_term_n_genes how many genes top genes from the PC should be inc. in the GO enrichment
#' @param n_PCA_top_genes how many genes from the PC should be plotted in a table
#' @param ntop how many highly variable genes should be inputted in the PCA
#' @param PC_cutoff minimal ammount of variance explained by PC for it to be plotted
#' @param gene_background optional list of background gene names
#' @param filename start of the output filename
#' @export
plotPCA_dds <- function (object, 
                              intgroup = c("condition1", "condition2"),
                              output_dir, 
                              run_go_enrichment = TRUE,
                              Genenames_keyType = 'SYMBOL', 
                              go_term_n_genes = 100, 
                              n_PCA_top_genes = 20, 
                              ntop = 1000,
                              PC_cutoff = 5, 
                              gene_background = 'not_provided',
                              filename = 'PCA_'){
  #object = deseq2 object to vizualize as a PCA plot
  #ingroup = specify the condition or condistions to vizualize, first one will be colour, second one will be shape
  #output_dir = location where the pdf files are saved
  #run_go_enrichment = TRUE/FALSE wether to run GO-term enrichment on the top 100 genes contributing to the PCs or not
  #Genenames_keyTYpe = enrichr 
  #go_term_n_genes = ammount of top PCA contributing genes GO-term-enrichment is ran upon.
  #n_PCA_top_genes = ammount of top PCA contributing genes plotted
  #ntop = ammount of variable genes that goes into the PCA plot.
  #PC_cutoff = percentage the PC needs to explain in order to be included in the graph
  #gene_background = optionally specify a list of gene names used as a background for GO-term enrichment. e.g. all genes measured with at least 10 counts across samples.
  
  rv <- matrixStats::rowVars(SummarizedExperiment::assay(object)) # calculate the variance for each gene
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))] # select the ntop genes by variance
  pca <- prcomp(t(SummarizedExperiment::assay(object)[select, ])) # perform a PCA on the data in assay(x) for the selected genes
  percentVar <- pca$sdev^2/sum(pca$sdev^2) # the contribution to the total variance for each component
  if (!all(intgroup %in% names(SummarizedExperiment::colData(object)))) {stop("the argument 'intgroup' should specify columns of colData(dds)")  }
  intgroup.df <- as.data.frame(SummarizedExperiment::colData(object)[, intgroup, drop = FALSE])
  if (length(colnames(intgroup.df)) == 1){
    print('only one intgroup column')
    group1 <- as.factor(SummarizedExperiment::colData(object)[[intgroup]])
    group2 <- as.factor(SummarizedExperiment::colData(object)[[intgroup]])}else {
    print('two intgroup columns')
    group1 <- as.factor(SummarizedExperiment::colData(object)[[intgroup[1]]])
    group2 <- as.factor(SummarizedExperiment::colData(object)[[intgroup[2]]]) }
  
  df_var <- as.data.frame(percentVar*100)
  df_var$PC <- paste('PC',row.names(df_var),sep="")
  df_var$PC <- factor(df_var$PC, levels = df_var$PC)
  colnames(df_var)<- c('percentage_variance','PC')
  var_plot <- ggplot2::ggplot(df_var, ggplot2::aes(x=PC, y = percentage_variance, group =1))+ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::geom_hline(yintercept = PC_cutoff)
  
  plot_list <- list()
  #find the ammount of PCs that explain more than 1% of the variance:
  sig_PCs <- nrow(df_var[df_var$percentage_variance >= PC_cutoff,])
  if (sig_PCs < 2){
    print(paste0('less than 2 PCs above variance cutoff of ',PC_cutoff,'%, taking 2 PCs as a backup'))
    sig_PCs <- 2
  }

  #make a combination matrix to plot all PCs to other PCs that explain more than PC_cutoff% 
  df <- expand.grid(c(1:sig_PCs),c(1:sig_PCs))
  df <- df[df$Var1 != df$Var2,]
  shuffle_PCs <- df$Var1 > df$Var2
  var1_shuffle <- df[shuffle_PCs,]$Var1
  var2_shuffle <- df[shuffle_PCs,]$Var2
  df[shuffle_PCs,]$Var1 <- var2_shuffle
  df[shuffle_PCs,]$Var2 <- var1_shuffle
  df <- unique(df)
  pca_rotation_df <- as.data.frame(pca$rotation)
  #extract PCA output to a df
  d <- as.data.frame(pca$x)
  d$group1 <- group1
  d$group2 <- group2
  d$name <-colnames(object)
  n_plot <- 0
  for (row in 1:nrow(df)) {
    n_plot <- n_plot + 1
    #generate the PC plot for this specific set of PCs
    PC1 <- paste('PC',df[row,]$Var1, sep = '')
    PC2 <- paste('PC',df[row,]$Var2, sep = '')
    print(PC1)
    print(PC2)
    varx <- as.integer(stringr::str_sub(PC1,-1, -1))
    vary <- as.integer(stringr::str_sub(PC2,-1, -1))
    plot_PCA <-ggplot2::ggplot(data = d, ggplot2::aes_string(x = PC1, y = PC2, color = group1, shape = group2)) +
      ggplot2::geom_point(size = 4) +
      #scale_shape_manual(values = c(3, 4, 7, 8, 15:18)) + #used shapes
      ggplot2::xlab(paste0(PC1, ": ", round(percentVar[varx] * 100), "% variance")) +
      ggplot2::ylab(paste0(PC2, ": ", round(percentVar[vary] * 100), "% variance")) +
      ggplot2::coord_fixed()
    #plot_list <- c(plot_list, as.grob(plot_PCA))
    #add the top contributing genes of each PC as a table next each plot:
    genes_contrib_df <- data.frame(PC1_high_genes = tail(row.names(pca_rotation_df[order(pca_rotation_df[,"PC1"]),]), n = n_PCA_top_genes),
                                   PC1_low_genes = head(row.names(pca_rotation_df[order(pca_rotation_df[,"PC1"]),]), n = n_PCA_top_genes),
                                   PC2_high_genes = tail(row.names(pca_rotation_df[order(pca_rotation_df[,"PC2"]),]), n = n_PCA_top_genes),
                                   PC2_low_genes = head(row.names(pca_rotation_df[order(pca_rotation_df[,"PC2"]),]), n = n_PCA_top_genes),
                                   stringsAsFactors=FALSE
    )
    colnames(genes_contrib_df) <-c(paste(PC1, "high"),paste(PC1, "low"),paste(PC2, "high"),paste(PC2, "low"))
      tt <- gridExtra::ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
    tbl <- gridExtra::tableGrob(genes_contrib_df, rows=NULL, theme=tt)
    #run GO-term enrichment on the top 200 genes contributing to each PC
    output_filename = paste0(output_dir,"/",filename,PC1,"_",PC2,'.pdf')

    if (run_go_enrichment){
      print('running GO')
      if (length(gene_background)>1){
        print('using costum gene background list')
        background <- gene_background}else{
        print('using standard org.HS.eg.db background for GO enrichment')
        background <- NA
        }
      
      PC_enrichment_plot_list = NULL 
      PC_enrichment_plot_list = list()
      plot_number = 1
      for (PC in c(PC1, PC2)){
        PC_enriched = FALSE
        GO_plot <-  grid::rectGrob(gp=grid::gpar(col=NA))
        print(paste0('running GO ',PC))
        print(paste0(PC,' high'))
        PC_high_ego <- clusterProfiler::enrichGO(gene         = row.names(head(pca_rotation_df[order(pca_rotation_df[,PC], decreasing = T),], go_term_n_genes)),
                            universe = background,
                            OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                            keyType       = Genenames_keyType,
                            ont           = "ALL",
                            pAdjustMethod = "BH",
                            qvalueCutoff  = 0.01,
                            )
        print(paste0(PC,' low'))
        PC_low_ego <- clusterProfiler::enrichGO(gene         = row.names(head(pca_rotation_df[order(pca_rotation_df[,PC], decreasing = F),], go_term_n_genes)),
                                universe = background,
                                 OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                 keyType       = Genenames_keyType,
                                 ont           = "ALL",
                                 pAdjustMethod = "BH",
                                 qvalueCutoff  = 0.01)
      #combine the plots if enrichment in both directions
      if ((nrow(PC_high_ego@result) > 0)& (nrow(PC_low_ego@result) > 0) ){
        print(paste0('combining plots ',PC))
        df_PC1_high <- as.data.frame(PC_high_ego)
        if (length(rownames(df_PC1_high)) > 10){df_PC1_high <- df_PC1_high[1:10,]}
        df_PC1_high$direction <- 1
        df_PC1_low <- as.data.frame(PC_low_ego)
        if (length(rownames(df_PC1_low)) > 10){df_PC1_low <- df_PC1_low[1:10,]}
        df_PC1_low$direction <- -1
        df_PC1 <- rbind(df_PC1_high, df_PC1_low)
        df_PC1$directed_count = df_PC1$direction * df_PC1$Count
        PC_enriched = TRUE
      }
      if(!PC_enriched){
        if ((nrow(PC_high_ego@result) > 0) & ((nrow(PC_low_ego@result) == 0))){
          print(paste0('high genes enriched ',PC))
          PC_high_ego <- as.data.frame(PC_high_ego)
          if (length(rownames(PC_high_ego)) > 10){PC_high_ego <- PC_high_ego[1:10,]}
          PC_high_ego$direction <- 1
          df_PC1 <- PC_high_ego
          df_PC1$directed_count = df_PC1$direction * df_PC1$Count
          PC_enriched = TRUE
        }}
      if(!PC_enriched){
        if ((nrow(PC_high_ego@result) == 0)& ( nrow(PC_low_ego@result) > 0) & !PC_enriched){
          print(paste0('low genes enriched ',PC))
          PC_low_ego <- as.data.frame(PC_low_ego)
          if (length(rownames(PC_low_ego)) > 10){PC_low_ego <- PC_low_ego[1:10,]}
          PC_low_ego$direction <- -1
          df_PC1 <- PC_low_ego
          df_PC1$directed_count = df_PC1$direction * df_PC1$Count
          PC_enriched = TRUE
        }}
      if (!PC_enriched){
        print( 'no GO terms enriched, trying out enrichment of all PC genes')
      }else{df_PC1 <- tidyr::drop_na(df_PC1)}
      if (length(df_PC1) > 0){
        GO_plot <- ggplot2::ggplot(data=df_PC1, ggplot2::aes(x=reorder(Description,directed_count), y=directed_count, fill = qvalue)) +
          ggplot2::geom_bar(stat="identity") + ggplot2::scale_fill_gradient(trans = "log", low = "red", high = "blue",breaks = c(0,0.000000001,0.000001, 0.001, 0.05, 0.5)) +
          ggplot2::coord_flip() + ggplot2::labs(x="", y = "" )}
      
      PC_enrichment_plot_list[[plot_number]] <-GO_plot 
      plot_number <- plot_number + 1
      }
      # Plot chart and table into one object
      print('adding finalplot to plotlist')
      final_plot <- grid::rectGrob(gp=grid::gpar(col=NA))
      #pdf(filename,width=24,height=12,paper='special', onefile = TRUE) 
      final_plot <- gridExtra::arrangeGrob(plot_PCA, PC_enrichment_plot_list[[2]],  
                                             ggplotify::as.grob(PC_enrichment_plot_list[[1]]), 
                                             ggplotify::as.grob(tbl),
                                             ncol=2,
                                             as.table=TRUE,
                                             heights=c(3,3))
      #dev.off()
      ggplot2::ggsave(file=output_filename, final_plot, dpi = 300, width = 40, height = 20, units = 'cm' ) #saves g
    }else{    # Plot chart and table into one object
      final_plot <- grid::rectGrob(gp=grid::gpar(col=NA))
      #final_plot <- grid.arrange(plot_PCA,
      final_plot <- gridExtra::arrangeGrob(plot_PCA,
                                 tbl,
                                 ncol=2,
                                 as.table=TRUE,
                                 heights=c(3,1))}
    #table_plot <- grid.table(genes_contrib_df)
    ggplot2::ggsave(file=output_filename, final_plot, dpi = 300, width = 40, height = 20, units = 'cm' ) #saves g
  }
  return(sig_PCs)
  
}
