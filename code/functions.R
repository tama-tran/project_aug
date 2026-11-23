## Functions
### Violin plot  # plot_vln_qc(seurat_list_1)
###

plot_vln_qc <- function(seurat_list){
  
  options(repr.plot.width = 14, repr.plot.height = 8)
  
  for(i in seq_along(seurat_list)){
    
    print(
      
      VlnPlot(
        seurat_list[[i]],
        features = c(
          "nFeature_RNA",
          "nCount_RNA",
          "percent.mt"
        ),
        ncol = 3
      ) +
        plot_annotation(
          title = names(seurat_list)[i]
        )
      
    )
    
  }
  
}

### Scatter plot # plot_scatter_qc(seurat_list_1)
###
plot_scatter_qc <- function(seurat_list){
  
  options(repr.plot.width = 16, repr.plot.height = 8)
  
  for(i in seq_along(seurat_list)){
    
    p1 <- FeatureScatter(
      seurat_list[[i]],
      feature1 = "nCount_RNA",
      feature2 = "nFeature_RNA"
    )
    
    p2 <- FeatureScatter(
      seurat_list[[i]],
      feature1 = "nCount_RNA",
      feature2 = "percent.mt"
    )
    
    print(p1 + p2)
    
  }
  
}

### Summary # qc_summary_1 <- get_qc_summary(seurat_list_1)
###
get_qc_summary <- function(seurat_list){
  
  do.call(
    rbind,
    lapply(seurat_list, function(obj){
      
      data.frame(
        
        Sample = unique(obj$sample),
        
        Sample_ID = unique(obj$sample_id),
        
        Cells = ncol(obj),
        
        Median_nFeature = median(obj$nFeature_RNA),
        
        Median_nCount = median(obj$nCount_RNA),
        
        Median_percent_mt = median(obj$percent.mt)
        
      )
      
    })
  )
  
}

### Subset      # seurat_list_1 <- filter_qc(seurat_list_1, qc_threshold_1)
###
filter_qc <- function(seurat_list, qc_threshold){
  
  for(i in seq_along(seurat_list)){
    
    seurat_list[[i]] <- subset(
      
      seurat_list[[i]],
      
      subset =
        nFeature_RNA >= qc_threshold$min_feature[i] &
        nFeature_RNA <= qc_threshold$max_feature[i] &
        nCount_RNA   <= qc_threshold$max_count[i] &
        percent.mt   <= qc_threshold$max_mt[i]
      
    )
    
  }
  
  return(seurat_list)
  
}

### workflow
#qc_before <- get_qc_summary(seurat_list_1)

#plot_vln_qc(seurat_list_1)

#plot_scatter_qc(seurat_list_1)

#qc_threshold <- data.frame(
#  max_feature = c(3000, 5000, 5000, 4000, 5000, 5000, 5000, 3500, 4000),
#  max_count   = c(18000, 20000, 20000, 18000, 20000, 20000, 20000, 18000, 20000),
#  max_mt      = c(25, 25, 28, 25, 25, 25, 25, 20, 20), #(20, 20, 28, 20, 20, 20, 20, 20, 20)
#  min_feature = c(300, 300, 300, 300, 300, 300, 300, 300, 600)
#                )
#seurat_list_1 <- filter_qc( seurat_list_1, qc_threshold_1)

#qc_after <- get_qc_summary(seurat_list_1)


############ DOWNSTREAM_ANALYSIS
build_panglaodb_list <- function(db,
                                 organ_filter = c("Liver", "Immune system"),
                                 min_spec     = 0.0) {
  if (!is.null(organ_filter)) {
    db <- db %>% filter(organ %in% organ_filter)
  }
  if (min_spec > 0) {
    db <- db %>% filter(!is.na(spec), spec >= min_spec)
  }
  marker_list <- db %>%
    group_by(cell_type) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")
  setNames(marker_list$genes, marker_list$cell_type)
}

#
build_cellmarker_list <- function(db, tissue_filter = "liver") {
  if (!is.null(tissue_filter)) {
    db <- db %>%
      filter(grepl(tissue_filter, tissue, ignore.case = TRUE))
  }
  marker_list <- db %>%
    group_by(cell_type) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")
  setNames(marker_list$genes, marker_list$cell_type)
}

#### core overlap functions
overlap_score <- function(cluster_genes,
                          marker_db,
                          top_n = 30) {
  
  cluster_genes <- head(cluster_genes, top_n)
  n_query       <- length(cluster_genes)
  
  if (n_query == 0) return(data.frame())
  
  #Gán kết quả pipe vào biến results
  results <- lapply(names(marker_db), function(cell_type) {
    ref_genes <- unique(marker_db[[cell_type]])
    n_ref     <- length(ref_genes)
    overlap   <- intersect(cluster_genes, ref_genes)
    n_overlap <- length(overlap)
    
    if (n_overlap == 0) return(NULL)
    
    union_size <- length(union(cluster_genes, ref_genes))
    jaccard    <- round(n_overlap / union_size, 4)
    pct_ref    <- round(n_overlap / n_ref   * 100, 1)
    pct_query  <- round(n_overlap / n_query * 100, 1)
    
    data.frame(
      cell_type         = cell_type,
      n_overlap         = n_overlap,
      n_ref_markers     = n_ref,
      pct_ref_covered   = pct_ref,
      pct_query_covered = pct_query,
      jaccard_index     = jaccard,
      matched_genes     = paste(overlap, collapse = ", "),
      stringsAsFactors  = FALSE
    )
  }) %>% 
    bind_rows()
  
  #Kiểm tra điều kiện ở dòng riêng biệt (không dùng %>% nữa)
  if (nrow(results) == 0) {
    return(data.frame())
  }
  
  #Sắp xếp kết quả thu được
  results <- results %>% 
    arrange(desc(jaccard_index))
  
  return(results)
}

### Run Across All Clusters and Sources
run_overlap_all_clusters <- function(top_markers, clusters,
                                     marker_db, source_label,
                                     top_n = 30) {
  lapply(clusters, function(cl) {
    genes  <- top_markers %>% filter(cluster == cl) %>% pull(gene)
    result <- overlap_score(genes, marker_db, top_n = top_n)
    if (is.null(result) || nrow(result) == 0) {
      return(data.frame(
        cell_type         = "Unknown",
        n_overlap         = 0,
        n_ref_markers     = 0,
        pct_ref_covered   = 0,
        pct_query_covered = 0,
        jaccard_index     = 0,
        matched_genes     = "",
        cluster           = cl,
        source            = source_label,
        stringsAsFactors  = FALSE))
    }
    result$cluster <- cl
    result$source  <- source_label
    result
  }) %>% bind_rows()
}


#####################
run_enrichment <- function(marker_df, type = "GO", ontology = "BP", universe_genes = NULL) {
  
  #Chuyển đổi SYMBOL sang ENTREZID cho các gene marker
  message("--> Converting gene SYMBOLs to ENTREZID...")
  gene.df <- bitr(
    unique(marker_df$gene),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
    drop = TRUE
  )
  
  # Nếu dùng universe (nền gene background), chuyển đổi cả universe
  universe_entrez <- NULL
  if (!is.null(universe_genes)) {
    univ.df <- bitr(universe_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
    universe_entrez <- univ.df$ENTREZID
  }
  
  # Map ENTREZID vào bảng marker ban đầu
  marker.df <- marker_df %>%
    left_join(gene.df, by = c("gene" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID))
  
  # Split gene theo từng cụm (cluster) và đảm bảo các gene trong cụm là unique
  gene_list <- split(marker.df$ENTREZID, marker.df$cluster)
  gene_list <- lapply(gene_list, unique)
  
  res_list <- list()
  
  # Loop qua từng cụm tế bào
  for(cl in names(gene_list)){
    message(paste0("Processing Cluster: ", cl, " | Type: ", type, if(type=="GO") paste0("-", ontology) else ""))
    
    genes_to_test <- gene_list[[cl]]
    
    # Bỏ qua nếu cụm này có quá ít gene marker
    if(length(genes_to_test) < 5) {
      message(paste0("Skipping cluster ", cl, " due to too few genes (< 5)."))
      next
    }
    
    #GO / KEGG
    if (type == "GO") {
      ego <- enrichGO(
        gene          = genes_to_test,
        universe      = universe_entrez,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE # Biến đổi ID thành tên Gene Symbol trực quan
      )
    } else if (type == "KEGG") {
      ego <- enrichKEGG(
        gene          = genes_to_test,
        universe      = universe_entrez,
        organism      = "hsa", # "hsa" là mã cho Homo sapiens (người)
        keyType       = "ncbi-geneid",  # "ncbi-geneid" là đại diện cho ENTREZID trong KEGG
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2
      )
      # Vì enrichKEGG không có tham số readable = TRUE trực tiếp, ta dùng hàm setReadable
      if(!is.null(ego) && nrow(as.data.frame(ego)) > 0){
        ego <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      }
    }

    if(!is.null(ego) && nrow(as.data.frame(ego)) > 0){
      df <- as.data.frame(ego)
      df$cluster <- cl
      res_list[[cl]] <- df
    }
  }

  if(length(res_list) > 0){
    return(do.call(rbind, res_list))
  } else {
    message("No enrichment pathways found for any cluster under the thresholds!")
    return(data.frame())
  }
} 