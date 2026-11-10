#########################################################
#########################################################
# Useful functions and bits
#########################################################
#########################################################

split_taxonomy <- function(counts_table) {
  tax_split <- strsplit(rownames(counts_table), ";")
  
  tax_levels <- data.frame(
    kingdom = gsub("k__", "", sapply(tax_split, `[`, 1)),
    phylum = gsub("p__", "", sapply(tax_split, `[`, 2)),
    class = gsub("c__", "", sapply(tax_split, `[`, 3)),
    order = gsub("o__", "", sapply(tax_split, `[`, 4)),
    family = gsub("f__", "", sapply(tax_split, `[`, 5)),
    genus = gsub("g__", "", sapply(tax_split, `[`, 6)),
    species = gsub("s__", "", sapply(tax_split, `[`, 7))
  )
  
  species_names <- tax_levels$species
  counts_data <- as.data.frame(lapply(counts_table, as.numeric))
  
  process_level <- function(data, group) {
    result <- aggregate(data, by = list(group), sum, na.rm = TRUE)
    rownames(result) <- result[,1]
    result[,-1]
  }
  
  list(
    kingdom = process_level(counts_data, tax_levels$kingdom),
    phylum = process_level(counts_data, tax_levels$phylum),
    class = process_level(counts_data, tax_levels$class),
    order = process_level(counts_data, tax_levels$order),
    family = process_level(counts_data, tax_levels$family),
    genus = process_level(counts_data, tax_levels$genus),
    species = process_level(counts_data, species_names)
  )
}

#########################################################
#########################################################
#########################################################

# Original filter_feature function
filter_feature <- function(count_df, tax_df = NULL,
                           filter_method = 'abs_count',
                           cutoff = 1, prev_cutoff = 2,
                           data_type = 'asv') {
  if(!data_type %in% c('asv','metagenomics')) {
    stop("data_type must be 'asv' or 'metagenomics'")
  }
  if(class(count_df) != 'data.frame') {
    stop("count_df must be a dataframe")
  }
  if(!filter_method %in% c('abs_count','percent_sample','percent_dataset')) {
    stop("filter_method must be 'abs_count','percent_sample', or 'percent_dataset'")
  }
  if(prev_cutoff > ncol(count_df)) {
    stop("prev_cutoff cannot exceed the number of samples")
  }
  if(filter_method == 'abs_count') {
    if(cutoff %% 1 != 0) {
      stop("cutoff must be an integer when filter_method = 'abs_count'")
    }
    if(cutoff > max(count_df)) {
      stop("cutoff cannot be greater than the maximum number of reads when filter_method = 'abs_count'")
    }
  } else {
    if(cutoff < 0 | cutoff > 100) {
      stop("cutoff must be between 0-100 when filter_method is 'percent_sample' or 'percent_dataset'")
    }
  }
  if(data_type == 'asv') {
    if(is.null(tax_df)) {
      stop("tax_df is required when data_type = 'asv'")
    }
    if(!identical(rownames(count_df), tax_df$featureID)) {
      stop("Order of features in count_df and tax_df must be identical")
    }
  }
  relab_by_sample <- apply(count_df, 2, function(x) x / sum(x))
  relab_by_data <- count_df / sum(count_df)
  prev_data <- rowSums(count_df != 0)
  filtered_count <- switch(
    filter_method,
    abs_count = count_df[rowSums(count_df >= cutoff) >= prev_cutoff,],
    percent_sample = count_df[rowSums(relab_by_sample >= cutoff) >= prev_cutoff,],
    percent_dataset = count_df[rowSums(relab_by_data >= cutoff) >= prev_cutoff,]
  )
  if(is.null(filtered_count)) {
    stop("No features remain after filtering. Please select different filter cutoff.")
  }
  sampleID <- colnames(filtered_count)
  empty_sample_ind <- which(colSums(filtered_count)==0)
  if(length(empty_sample_ind) > 0) {
    filtered_count <- filtered_count[, sampleID[-empty_sample_ind]]
  }
  featID <- rownames(filtered_count)
  empty_feat_ind <- which(rowSums(filtered_count)==0)
  if(length(empty_feat_ind) > 0) {
    filtered_count <- filtered_count[featID[-empty_feat_ind],]
  }
  feature_keep <- rownames(filtered_count)
  feature_remove <- setdiff(rownames(count_df), rownames(filtered_count))
  if(data_type == 'asv') {
    tax_filtered <- tax_df %>%
      filter(featureID %in% feature_keep)
  } else {
    tax_filtered <- NULL
  }
  res <- list(
    'filtered' = filtered_count,
    'feature_remove' = feature_remove,
    'feature_keep' = feature_keep
  )
  if(data_type == 'asv') {
    res$taxonomy <- tax_filtered
  }
  return(res)
}

#########################################################
#########################################################
#########################################################

# Function to run DESeq2 analysis
run_deseq <- function(txi, metadata, design_formula) {
  dds <- DESeqDataSetFromTximport(txi, metadata, design = design_formula)
  dds <- dds[rowSums(counts(dds) >= 10) >= 15, ]
  dds <- DESeq(dds, test = "Wald")
  res <- as.data.frame(results(dds))
  return(list(dds = dds, res = res))
}

#########################################################
#########################################################
#########################################################

# Function to annotate DESeq2 results with gene names
process_deseq_results <- function(res, mart = NULL, offline = FALSE) {
  gene_ids <- gsub("(\\.).*", "", rownames(res))
  
  if (offline) {
    # Use org.Mm.eg.db for offline annotation
    cat("Using offline annotation with org.Mm.eg.db\n")
    
    # More efficient mapping with direct key lookup
    ensembl_to_symbol <- suppressWarnings(AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = gene_ids,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    ))
    
    # Create a more efficient data frame
    gene_map <- data.frame(
      gene_id = names(ensembl_to_symbol),
      gene_name = unname(ensembl_to_symbol),
      stringsAsFactors = FALSE
    )
  } else {
    # Use biomaRt for online annotation
    cat("Using online annotation with biomaRt\n")
    symbols <- getBM(
      attributes = c('external_gene_name', 'ensembl_gene_id'),
      filters = 'ensembl_gene_id',
      values = gene_ids,
      mart = mart
    )
    gene_map <- data.frame(gene_id = symbols$ensembl_gene_id, gene_name = symbols$external_gene_name)
  }
  # Create a named vector for faster lookup
  name_lookup <- setNames(gene_map$gene_name, gene_map$gene_id)
  
  # Add gene names to results more efficiently
  res$gene_id <- gene_ids
  res$gene_name <- name_lookup[res$gene_id]
  
  duplicate_count <- sum(duplicated(res$gene_name[!is.na(res$gene_name)]))
  cat("Number of duplicate gene names removed:", duplicate_count, "\n")
  
  res <- res[!duplicated(res$gene_name) & !is.na(res$gene_name), ]
  rownames(res) <- res$gene_name
  res <- subset(res, !grepl("^Gm", gene_name))
  
  symbol_to_ensembl <- setNames(gene_map$gene_id, gene_map$gene_name)
  
  return(list(res = res, symbol_to_ensembl = symbol_to_ensembl, duplicate_count = duplicate_count))
}

#########################################################
#########################################################
#########################################################

generate_volcano_plot <- function(res, title) {
  # Create a new column to categorize genes
  res$category <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, 
                         "pvalue + log2FC", "NS")
  
  # Convert p-values to -log10 scale
  res$neg_log10_padj <- -log10(res$padj)
  
  # Calculate cutoff line for p-value
  p_cutoff <- -log10(0.05)
  
  # Create the plot
  ggplot(res, aes(x = log2FoldChange, y = neg_log10_padj, color = category)) +
    # Add threshold lines
    geom_hline(yintercept = p_cutoff, linetype = "dashed", color = "black", size = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5) +
    # Add points
    geom_point(size = 1.0, alpha = 1) +
    # Set color scheme
    scale_color_manual(values = c("NS" = "grey", "pvalue + log2FC" = "red")) +
    # Labels
    labs(
      x = expression(Log[2]~"fold change"),
      y = expression((-log[10])~"adjusted p-value"),
      title = title
    ) +
    # Theme adjustments to match ComplexHeatmap style
    theme_classic() +
    theme(
      text = element_text(size = 20),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      legend.position = "right",
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black", size = 0.3),
      legend.key.size = unit(1, "cm"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.title = element_text(hjust = 0.5)
    )
}

#########################################################
#########################################################
#########################################################

generate_ridge_plot <- function(fgsea_res, title) {
  top_terms <- fgsea_res[order(fgsea_res$padj), ]
  if (nrow(top_terms) > 10) {
    top_terms <- top_terms[1:10, ]
  }
  
  top_terms$pathway <- gsub("GOBP_", "", top_terms$pathway)
  top_terms$pathway <- gsub("_", " ", top_terms$pathway)
  
  clean_title <- gsub("_", " ", title)
  
  plot <- ggplot(top_terms, aes(x = NES, y = reorder(pathway, NES))) +
    geom_bar(stat = "identity", aes(fill = padj), width = 0.9) +
    scale_fill_gradient(low = "red", high = "blue") +
    theme_minimal() +
    labs(title = clean_title, x = "Normalized Enrichment Score", y = "Pathway") +
    theme(
      text = element_text(size = 10),
      axis.text.y = element_text(size = 9, lineheight = 0.9),
      axis.text.x = element_text(size = 9),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 11),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      legend.key.size = unit(0.8, 'lines'),
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(5, 10, 5, 20),
      panel.grid.major.y = element_blank(),
      aspect.ratio = 0.7
    ) +
    scale_y_discrete(labels = function(x) {
      sapply(strwrap(x, width = 30, simplify = FALSE),
             paste, collapse = "\n")
    })
  
  return(list(plot = plot, top_terms = top_terms))
}

#########################################################
#########################################################
#########################################################

# Running DESeq
run_deseq_analysis <- function(comparisons, mart = NULL, offline = FALSE) {
  deseq_results <- list()
  
  for (comp in comparisons) {
    cat("Running DESeq2 for comparison:", comp$name, "\n")
    
    # Extract metadata and design formula from comparison
    deseq_out <- run_deseq(comp$txi, comp$metadata, comp$design_formula)
    
    # Process DESeq2 results
    processed_res <- process_deseq_results(deseq_out$res, mart, offline)
    annotated_res <- processed_res$res
    symbol_to_ensembl <- processed_res$symbol_to_ensembl
    
    # Generate volcano plot
    volcano_plot <- generate_volcano_plot(annotated_res, comp$title)
    
    # Store DESeq2 results
    deseq_results[[comp$name]] <- list(
      deseq = annotated_res,
      dds = deseq_out$dds,
      symbol_to_ensembl = symbol_to_ensembl,
      volcano_plot = volcano_plot,
      duplicate_count = processed_res$duplicate_count
    )
  }
  
  return(deseq_results)
}

#########################################################
#########################################################
#########################################################

run_fgsea_analysis <- function(deseq_results, gmt_path, 
                               pval_cutoff = 0.05) {
  # Load pathways just once
  pathways <- gmtPathways(gmt_path)
  gsea_results <- list()
  
  for (comp_name in names(deseq_results)) {
    cat("Running fgsea for comparison:", comp_name, "\n")
    annotated_res <- deseq_results[[comp_name]]$deseq
    
    # Extract gene fold changes more efficiently
    genes_fc <- annotated_res$log2FoldChange
    names(genes_fc) <- rownames(annotated_res)
    
    # Sort genes by fold change for faster fgsea processing
    genes_fc <- sort(genes_fc, decreasing = TRUE)
    
    
    fgsea_res <- fgsea(pathways, genes_fc)
    
    
    significant_pathways <- fgsea_res[fgsea_res$padj < pval_cutoff, ]
    
    if (nrow(significant_pathways) > 0) {
      ridge_plot_list <- generate_ridge_plot(significant_pathways, paste0(comp_name, " Pathways"))
      ridge_plot <- ridge_plot_list$plot
      top_terms <- ridge_plot_list$top_terms
    } else {
      top_pathways <- fgsea_res[order(-abs(fgsea_res$NES)), ][1:10, ]
      ridge_plot_list <- generate_ridge_plot(top_pathways, paste0(comp_name, " Top Pathways (Not Significant)"))
      ridge_plot <- ridge_plot_list$plot
      top_terms <- ridge_plot_list$top_terms
      cat("No significant pathways found for", comp_name, "at padj <", pval_cutoff, "\n")
    }
    
    gsea_results[[comp_name]] <- list(
      gsea = fgsea_res,
      significant_pathways = significant_pathways,
      ridge_plot = ridge_plot,
      top_terms = top_terms,
      leading_edge_genes = if(nrow(significant_pathways) > 0) { 
        setNames(significant_pathways$leadingEdge, significant_pathways$pathway)
      } else {
        NULL
      }
    )
  }
  
  return(gsea_results)
}

#########################################################
#########################################################
#########################################################

run_full_analysis <- function(comparisons, mart = NULL, offline = FALSE, 
                              gmt_path = "/well/powrie/shared/mm_msigdb/m5.go.v2023.1.Mm.symbols.gmt.txt",
                              go_category = "BP") {
  # Step 1: Run DESeq2 analysis
  cat("Running DESeq2 analysis...\n")
  deseq_results <- run_deseq_analysis(comparisons, mart, offline)
  
  # Step 2: Filter GMT file for specific GO category and save to temp file
  cat(paste0("Filtering GMT for GO", go_category, " terms...\n"))
  temp_gmt_path <- tempfile(pattern = paste0("GO", go_category, "_filtered_"), fileext = ".gmt")
  
  # Read the GMT file
  gmt_lines <- readLines(gmt_path)
  
  # Filter lines that start with the specified GO category
  category_prefix <- paste0("GO", go_category, "_")
  filtered_lines <- gmt_lines[grepl(paste0("^", category_prefix), gmt_lines)]
  
  # Write filtered lines to temp file
  writeLines(filtered_lines, temp_gmt_path)
  
  # Step 3: Run fgsea analysis with filtered pathways
  cat("Running fgsea analysis...\n")
  gsea_results <- run_fgsea_analysis(deseq_results, gmt_path = temp_gmt_path)
  
  # Clean up temp file
  file.remove(temp_gmt_path)
  
  # Combine results
  results <- list()
  for (comp_name in names(deseq_results)) {
    results[[comp_name]] <- c(
      deseq_results[[comp_name]],
      gsea_results[[comp_name]]
    )
  }
  
  return(results)
}

#########################################################
#########################################################
#########################################################

# Function to properly create txi objects
create_txi_subset <- function(txi, samples) {
  new_txi <- list(
    abundance = txi$abundance[, samples],
    counts = txi$counts[, samples],
    length = txi$length[, samples],
    countsFromAbundance = txi$countsFromAbundance
  )
  # Preserve any other attributes
  attributes(new_txi) <- attributes(txi)
  return(new_txi)
}

#########################################################
#########################################################
#########################################################

perform_kruskal_dunn_tests <- function(data) {
  numeric_columns <- sapply(data, is.numeric)
  
  kruskal_results <- lapply(names(data)[numeric_columns], function(col) {
    test_result <- kruskal.test(data[[col]] ~ data$status)
    data.frame(
      Variable = col,
      Chi_Squared = test_result$statistic,
      Df = test_result$parameter,
      P_Value = test_result$p.value
    )
  })
  
  kruskal_df <- do.call(rbind, kruskal_results)
  kruskal_df$P_Adjusted <- p.adjust(kruskal_df$P_Value, method = "BH")
  
  post_hoc_dunn <- function(column_name) {
    dunn_test <- dunnTest(data[[column_name]] ~ data$status, method = "bh")
    dunn_df <- dunn_test$res
    dunn_df$Variable <- column_name
    
    # Determine which group has the higher median
    groups <- levels(data$status)
    medians <- aggregate(data[[column_name]], by = list(data$status), median)
    names(medians) <- c("Group", "Median")
    
    dunn_df$Higher_Group <- apply(dunn_df, 1, function(row) {
      group1 <- unlist(strsplit(as.character(row["Comparison"]), " - "))[1]
      group2 <- unlist(strsplit(as.character(row["Comparison"]), " - "))[2]
      median1 <- medians$Median[medians$Group == group1]
      median2 <- medians$Median[medians$Group == group2]
      if (median1 > median2) {
        return(group1)
      } else {
        return(group2)
      }
    })
    
    return(dunn_df)
  }
  
  significant_vars <- kruskal_df$Variable[kruskal_df$P_Adjusted < 0.05]
  post_hoc_results <- lapply(significant_vars, post_hoc_dunn)
  post_hoc_df <- do.call(rbind, post_hoc_results)
  
  list(kruskal_df = kruskal_df, post_hoc_df = post_hoc_df)
}

#########################################################
#########################################################
#########################################################

run_adonis2_tests <- function(normed_dat, metadata, factors, method = "euclidean", seed = 123) {
  results <- list()
  
  set.seed(seed)
  results[["status"]] <- adonis2(t(normed_dat) ~ status, data = metadata, method = method)
  
  for (factor in factors) {
    set.seed(seed)
    formula_str <- paste0("t(normed_dat) ~ status + `", factor, "`")
    results[[factor]] <- adonis2(as.formula(formula_str), data = metadata, method = method)
  }
  
  return(results)
}

#########################################################
#########################################################
#########################################################

create_adonis_df <- function(adonis_results, factors) {
  data.frame(
    Factor = c("status", factors),
    Df = c(adonis_results[["status"]]$Df[1], 
           sapply(adonis_results[factors], function(x) x$Df[2])),
    SumOfSqs = c(adonis_results[["status"]]$SumOfSqs[1], 
                 sapply(adonis_results[factors], function(x) x$SumOfSqs[2])),
    R2 = c(adonis_results[["status"]]$R2[1], 
           sapply(adonis_results[factors], function(x) x$R2[2])),
    F = c(adonis_results[["status"]]$F[1], 
          sapply(adonis_results[factors], function(x) x$F[2])),
    P_value = c(adonis_results[["status"]]$`Pr(>F)`[1], 
                sapply(adonis_results[factors], function(x) x$`Pr(>F)`[2]))
  )
}

#########################################################
#########################################################
#########################################################

convert_to_counts <- function(transcript_mapping_file, input_directory) {
  tx2gene <- read.csv(transcript_mapping_file, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  currwd <- getwd()
  setwd(input_directory)
  files <- list.files()[grep("*abundance.tsv", list.files())]
  names(files) <- gsub("_abundance.tsv", "", files)
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance="lengthScaledTPM")
  counts <- txi$counts
  setwd(currwd)
  row.names(counts) <- gsub("\\..*", "", row.names(counts))
  write.table(counts, "counts.tsv", sep = '\t')
  return(counts)
}

#########################################################
#########################################################
#########################################################

aggregate_counts_by_gene_name <- function(counts, ensembl_dataset) {
  ensembl <- useMart('ensembl', dataset = ensembl_dataset)
  attributes = c("ensembl_gene_id", "external_gene_name")
  gene_info = getBM(attributes = attributes,
                    filters = 'ensembl_gene_id', 
                    values = row.names(counts),
                    mart = ensembl,
                    mirror = 'useast')
  row.names(gene_info) <- gene_info$ensembl_gene_id
  counts <- counts[gene_info$ensembl_gene_id,]
  counts <- aggregate(counts, by=list(gene_info$external_gene_name), FUN='sum')
  row.names(counts) <- counts$Group.1
  counts <- counts[,2:ncol(counts)]
  counts <- as.matrix(counts)
  counts <- counts[-1,]
  return(counts)
}

#########################################################
#########################################################
#########################################################

local_aggregate_counts_by_gene_name <- function(counts, local_dataset) {
  row.names(counts) <- gsub("\\..*", "", row.names(counts))
  gene_info <- read.delim('/well/powrie/shared/music_examples/mm39_ensembl107_gene_id2gene_name.tsv', sep = '\t', header = T)
  row.names(gene_info) <- gene_info$ensemble_gene_id
  
  genes <- intersect(rownames(counts), gene_info$ensemble_gene_id)
  counts <- counts[genes,]
  gene_info <- gene_info[genes,]
  
  counts <- aggregate(counts, by=list(gene_info$external_ensemble_gene_id), FUN='sum')
  row.names(counts) <- counts$Group.1
  counts <- counts[,2:ncol(counts)]
  counts <- as.matrix(counts)
  counts <- counts[-1,]
  return(counts)
}

#########################################################
#########################################################
#########################################################

bulkdata <- function(counts){
  counts <- as.matrix(counts)
  bulkexp_set <- ExpressionSet(assayData = counts)
  bulk.mtx = exprs(bulkexp_set)
  return(bulk.mtx)
}

#########################################################
#########################################################
#########################################################

readSCE <- function(file){
  sce <- readRDS(file)
  return(sce)
}

#########################################################
#########################################################
#########################################################

get_top_loading_genes <- function(pca_result, pc_number = 1, n_top_genes = 15) {
  loadings <- pca_result$rotation[, pc_number]
  abs_loadings <- as.data.frame(abs(loadings))
  colnames(abs_loadings) <- "loading"
  abs_loadings$ensembl <- rownames(pca_result$rotation)
  
  gene_symbols <- mapIds(
    org.Mm.eg.db,
    keys = abs_loadings$ensembl,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  abs_loadings$gene_symbol <- gene_symbols
  abs_loadings <- abs_loadings[!is.na(abs_loadings$gene_symbol), ]
  abs_loadings <- abs_loadings[order(abs_loadings$loading, decreasing = TRUE), ]
  
  top_genes <- abs_loadings$ensembl[1:min(n_top_genes, nrow(abs_loadings))]
  
  return(list(
    top_genes = top_genes,
    abs_loadings = abs_loadings
  ))
}

#########################################################
#########################################################
#########################################################

get_top_loading_genes_pos_neg <- function(pca_result, pc_number = 1, n_top_genes = 15) {
  loadings <- pca_result$rotation[, pc_number]
  loadings_df <- data.frame(
    loading = loadings,
    ensembl = rownames(pca_result$rotation)
  )
  
  # Get gene symbols
  gene_symbols <- mapIds(
    org.Mm.eg.db,
    keys = loadings_df$ensembl,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  loadings_df$gene_symbol <- gene_symbols
  loadings_df <- loadings_df[!is.na(loadings_df$gene_symbol), ]
  
  # Split into positive and negative loadings
  positive_loadings <- loadings_df[loadings_df$loading > 0, ]
  negative_loadings <- loadings_df[loadings_df$loading < 0, ]
  
  # Order positive loadings (highest first)
  positive_loadings <- positive_loadings[order(positive_loadings$loading, decreasing = TRUE), ]
  
  # Order negative loadings (most negative first)
  negative_loadings <- negative_loadings[order(negative_loadings$loading, decreasing = FALSE), ]
  
  # Get top genes
  n_pos <- min(n_top_genes, nrow(positive_loadings))
  n_neg <- min(n_top_genes, nrow(negative_loadings))
  
  top_positive_genes <- if(n_pos > 0) positive_loadings$ensembl[1:n_pos] else character(0)
  top_negative_genes <- if(n_neg > 0) negative_loadings$ensembl[1:n_neg] else character(0)
  
  return(list(
    top_positive_genes = top_positive_genes,
    top_negative_genes = top_negative_genes,
    positive_loadings = positive_loadings[1:n_pos, ],
    negative_loadings = negative_loadings[1:n_neg, ]
  ))
}

#########################################################
#########################################################
#########################################################


plot_gene_expression_heatmap <- function(expression_data, top_genes, abs_loadings, 
                                         metadata, sample_id_col = "sampleID", 
                                         group_col = "status", title = "Gene Expression Heatmap") {
  
  top_gene_expr <- expression_data[top_genes, , drop = FALSE]
  gene_expr_long <- reshape2::melt(top_gene_expr)
  colnames(gene_expr_long) <- c("ENSEMBL", "Sample", "Expression")
  
  gene_expr_long$Gene <- abs_loadings$gene_symbol[match(gene_expr_long$ENSEMBL, 
                                                        abs_loadings$ensembl)]
  
  sample_data <- metadata[match(gene_expr_long$Sample, metadata[[sample_id_col]]), ]
  gene_expr_long[[group_col]] <- sample_data[[group_col]]
  
  heatmap_plot <- ggplot(gene_expr_long, aes(x = Sample, y = Gene, fill = Expression)) +
    geom_tile() +
    scale_fill_viridis_c(option = 'D') +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(title = title, x = "", y = "") +
    facet_grid(. ~ get(group_col), scales = "free_x", space = "free_x")
  
  return(heatmap_plot)
}

#########################################################
#########################################################
#########################################################

plot_pc_gene_heatmaps <- function(pca_result, expression_data, metadata, 
                                  n_top_genes = 15, sample_id_col = "sampleID", 
                                  group_col = "status") {
  
  pc1_result <- get_top_loading_genes(pca_result, pc_number = 1, n_top_genes = n_top_genes)
  
  pc2_result <- get_top_loading_genes(pca_result, pc_number = 2, n_top_genes = n_top_genes)
  
  pc1_heatmap <- plot_gene_expression_heatmap(
    expression_data = expression_data,
    top_genes = pc1_result$top_genes,
    abs_loadings = pc1_result$abs_loadings,
    metadata = metadata,
    sample_id_col = sample_id_col,
    group_col = group_col,
    title = "Expression of top PC1 genes"
  )
  
  pc2_heatmap <- plot_gene_expression_heatmap(
    expression_data = expression_data,
    top_genes = pc2_result$top_genes,
    abs_loadings = pc2_result$abs_loadings,
    metadata = metadata,
    sample_id_col = sample_id_col,
    group_col = group_col,
    title = "Expression of top PC2 genes"
  )
  
  return(list(
    pc1_heatmap = pc1_heatmap,
    pc2_heatmap = pc2_heatmap
  ))
}


#########################################################
#########################################################
#########################################################

get_combined_pathways <- function(pathways_df, n = 10, tissue_name) {
  df <- data.table::copy(pathways_df)
  
  spf_enriched <- df[NES > 0][order(-NES)][1:n]
  spf_enriched$Enriched_in <- "SPF"
  
  gf_enriched <- df[NES < 0][order(NES)][1:n]
  gf_enriched$Enriched_in <- "GF"
  
  combined <- rbindlist(list(spf_enriched, gf_enriched))
  combined$Tissue <- tissue_name
  
  combined$pathway <- gsub("GOBP_", "", combined$pathway)
  combined$pathway <- gsub("_", " ", combined$pathway)
  
  return(combined)
}

#########################################################
#########################################################
#########################################################

compute_recovery <- function(up_genes, down_genes, mm12_up_genes, mm12_down_genes, tissue) {
  # For upregulated genes
  up_not_recovered <- intersect(up_genes$gene_name, mm12_up_genes$gene_name)
  up_recovered <- setdiff(up_genes$gene_name, mm12_up_genes$gene_name)
  
  n_not_recovered_up <- length(up_not_recovered)
  n_recovered_up <- length(up_genes$gene_name) - n_not_recovered_up
  
  p_recovered_up <- (n_recovered_up/length(up_genes$gene_name))*100
  p_not_recovered_up <- (n_not_recovered_up/length(up_genes$gene_name))*100
  
  recovery_up <- data.frame(proportion = c(p_recovered_up, p_not_recovered_up),
                            recovery_status = c("Recovered", "Not recovered"),
                            tissue = rep(tissue, 2),
                            group = rep(paste(tissue, "up"), 2))
  
  # For downregulated genes
  down_not_recovered <- intersect(down_genes$gene_name, mm12_down_genes$gene_name)
  down_recovered <- setdiff(down_genes$gene_name, mm12_down_genes$gene_name)
  
  n_not_recovered_down <- length(down_not_recovered)
  n_recovered_down <- length(down_genes$gene_name) - n_not_recovered_down
  
  p_recovered_down <- (n_recovered_down/length(down_genes$gene_name))*100
  p_not_recovered_down <- (n_not_recovered_down/length(down_genes$gene_name))*100
  
  recovery_down <- data.frame(proportion = c(p_recovered_down, p_not_recovered_down),
                              recovery_status = c("Recovered", "Not recovered"),
                              tissue = rep(tissue, 2),
                              group = rep(paste(tissue, "down"), 2))
  
  list(up = recovery_up, down = recovery_down)
}

#########################################################
#########################################################
#########################################################

get_recovery_genes <- function(up_genes, down_genes, mm12_up_genes, mm12_down_genes, tissue) {
  # For upregulated genes
  up_not_recovered <- intersect(up_genes$gene_name, mm12_up_genes$gene_name)
  up_recovered <- setdiff(up_genes$gene_name, mm12_up_genes$gene_name)
  
  # For downregulated genes
  down_not_recovered <- intersect(down_genes$gene_name, mm12_down_genes$gene_name)
  down_recovered <- setdiff(down_genes$gene_name, mm12_down_genes$gene_name)
  
  return(list(up_recovered = up_recovered,
              down_recovered = down_recovered,
              up_not_recovered = up_not_recovered,
              down_not_recovered = down_not_recovered))
}


#########################################################
#########################################################
#########################################################

extract_wgcna_stats <- function(expr_data, module_colors, tissue_name) {
  # Get gene names from expression data
  gene_names <- colnames(expr_data)
  
  # Make sure module_colors has the same length as number of genes
  if(length(module_colors) != length(gene_names)) {
    stop("Length of module_colors doesn't match number of genes")
  }
  
  # Create named vector
  names(module_colors) <- gene_names
  
  # Get unique modules (excluding grey)
  unique_modules <- unique(module_colors[module_colors != "grey"])
  
  module_stats <- data.frame()
  
  for(i in seq_along(unique_modules)) {
    module <- unique_modules[i]
    module_genes <- gene_names[module_colors == module]
    module_expr <- expr_data[, module_genes, drop = FALSE]
    
    gene_sds <- apply(module_expr, 2, sd)
    all_expr_values <- as.vector(as.matrix(module_expr))
    
    stats <- data.frame(
      Module = paste0(tissue_name, i, " - ", module),
      Number_of_Genes = length(module_genes),
      Mean_Gene_SD = round(mean(gene_sds), 3),
      Median_Gene_SD = round(median(gene_sds), 3),
      Mean_Expression = round(mean(all_expr_values), 3),
      Expression_SD = round(sd(all_expr_values), 3)
    )
    
    module_stats <- rbind(module_stats, stats)
  }
  
  module_stats <- module_stats[order(-module_stats$Number_of_Genes), ]
  return(module_stats)
}


#########################################################
#########################################################
#########################################################

analyzeModuleData <- function(mergedColors, count_data, tissue_name = "Tissue") {
  # Define helper function within main function
  convertEnsemblToSymbol = function(ensembl_ids) {
    ensembl_ids = gsub("\\.[0-9]+$", "", ensembl_ids)
    
    gene_symbols = mapIds(org.Mm.eg.db,
                          keys = ensembl_ids,
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")
    
    gene_symbols[is.na(gene_symbols)] = ensembl_ids[is.na(gene_symbols)]
    
    return(gene_symbols)
  }
  
  # Load all required libraries
  library(circlize)
  library(ComplexHeatmap)
  library(reshape2)
  library(pheatmap)
  
  # Store all plots in a list
  all_plots <- list()
  
  modules = unique(mergedColors)
  all_go_results = list()
  
  gene_names = colnames(count_data)
  
  if(length(gene_names) != length(mergedColors)) {
    stop(paste("Gene names from", tissue_name, "count data don't match length of mergedColors"))
  }
  
  if(is.null(names(mergedColors))) {
    names(mergedColors) = gene_names
  }
  
  for(module in modules) {
    if(module == "grey") next
    
    module_indices = which(mergedColors == module)
    
    if(length(module_indices) == 0) {
      cat("No genes found in module", module, "\n")
      next
    }
    
    module_genes = gene_names[module_indices]
    
    cat("Processing module", module, "with", length(module_genes), "genes\n")
    
    go_result = try(enrichGO(gene = module_genes,
                             OrgDb = org.Mm.eg.db,
                             keyType = "SYMBOL",
                             ont = "BP",
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05,
                             readable = TRUE))
    
    if(!inherits(go_result, "try-error") && !is.null(go_result) && nrow(go_result@result) > 0) {
      # Store plots in the list
      all_plots[[paste0(module, "_dotplot")]] <- dotplot(go_result, showCategory = 15, 
                                                         title = paste("GO Terms -", module, "Module -", tissue_name))
      
      all_plots[[paste0(module, "_barplot")]] <- barplot(go_result, showCategory = 15, 
                                                         title = paste("GO Terms -", module, "Module -", tissue_name))
      
      all_go_results[[module]] = go_result
    } else {
      cat("No significant GO terms found for module", module, "\n")
    }
  }
  
  summary_df = data.frame(Module = character(),
                          GOterm = character(),
                          Pvalue = numeric(),
                          GeneRatio = character(),
                          Count = numeric(),
                          stringsAsFactors = FALSE)
  
  for(module in names(all_go_results)) {
    if(nrow(all_go_results[[module]]@result) > 0) {
      top_go = all_go_results[[module]]@result[1:min(5, nrow(all_go_results[[module]]@result)),]
      
      for(i in 1:nrow(top_go)) {
        summary_df = rbind(summary_df, data.frame(
          Module = module,
          GOterm = top_go$Description[i],
          Pvalue = top_go$pvalue[i],
          GeneRatio = top_go$GeneRatio[i],
          Count = top_go$Count[i],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  if(nrow(summary_df) > 0) {
    pvalue_matrix = acast(summary_df, GOterm ~ Module, value.var = "Pvalue", fill = 1)
    
    pvalue_matrix = -log10(pvalue_matrix)
    top_terms = names(sort(apply(pvalue_matrix, 1, max), decreasing = TRUE)[1:min(50, nrow(pvalue_matrix))])
    
    pvalue_matrix_filtered = pvalue_matrix[top_terms, ]
    
    # Store the pheatmap in the list
    all_plots[["go_term_heatmap"]] <- pheatmap(pvalue_matrix, 
                                               main = paste(tissue_name, "- Module vs GO Term Enrichment"),
                                               color = colorRampPalette(c("pink", "purple"))(50),
                                               fontsize = 15)
  }
  
  if(nrow(summary_df) > 0) {
    lowest_pvalue_terms <- aggregate(Pvalue ~ Module + GOterm, summary_df, min)
    lowest_pvalue_per_module <- aggregate(Pvalue ~ Module, lowest_pvalue_terms, min)
    
    module_terms <- merge(lowest_pvalue_per_module, lowest_pvalue_terms)
    
    all_modules <- unique(summary_df$Module)
    all_selected_terms <- unique(module_terms$GOterm)
    
    pvalue_matrix <- matrix(1, 
                            nrow = length(all_selected_terms), 
                            ncol = length(all_modules),
                            dimnames = list(all_selected_terms, all_modules))
    
    for(i in 1:nrow(module_terms)) {
      row <- module_terms$GOterm[i]
      col <- module_terms$Module[i]
      pvalue_matrix[row, col] <- module_terms$Pvalue[i]
    }
    
    pvalue_matrix <- -log10(pvalue_matrix)
    
    # Store the pheatmap in the list
    all_plots[["top_go_term_heatmap"]] <- pheatmap(pvalue_matrix, 
                                                   main = paste(tissue_name, "- Module vs Top GO Term Enrichment (Lowest p-value per module)"),
                                                   color = colorRampPalette(c("pink", "purple"))(50),
                                                   fontsize = 15)
  }
  
  plotModuleTopGenesWithGO = function(module) {
    cat("\n==== Processing module:", module, "====\n")
    
    module_indices = which(mergedColors == module)
    
    if(length(module_indices) == 0) {
      cat("No genes found in module", module, "\n")
      return(NULL)
    }
    
    module_genes = names(mergedColors)[module_indices]
    valid_genes = intersect(module_genes, colnames(count_data))
    
    if(length(valid_genes) == 0) {
      cat("None of the module genes found in expression data\n")
      return(NULL)
    }
    
    module_expr = t(count_data[, valid_genes, drop=FALSE])
    
    gene_vars = apply(module_expr, 1, var)
    valid_vars = which(!is.na(gene_vars) & gene_vars > 0)
    
    if(length(valid_vars) == 0) {
      cat("No genes with valid variance\n")
      return(NULL)
    }
    
    valid_genes = rownames(module_expr)[valid_vars]
    sorted_genes = valid_genes[order(gene_vars[valid_vars], decreasing = TRUE)]
    top_genes = sorted_genes[1:min(30, length(sorted_genes))]
    
    top_expr = module_expr[top_genes, , drop=FALSE]
    
    if(any(grepl("^ENS", rownames(top_expr)))) {
      gene_symbols = convertEnsemblToSymbol(rownames(top_expr))
      rownames(top_expr) = gene_symbols
    }
    
    top_expr_scaled = t(scale(t(as.matrix(top_expr))))
    
    go_annotation = "No significant GO terms"
    plot_title = paste("Module", module)
    
    if(exists("summary_df")) {
      module_go_terms = summary_df[summary_df$Module == module, "GOterm"]
      if(length(module_go_terms) > 0) {
        top_go_terms = module_go_terms[1:min(3, length(module_go_terms))]
        go_annotation = paste(module_go_terms[1:min(5, length(module_go_terms))], collapse = "\n")
        plot_title = paste(top_go_terms, collapse = "\n")
      }
    }
    
    col_fun = colorRamp2(
      c(min(top_expr_scaled, na.rm = TRUE), 0, max(top_expr_scaled, na.rm = TRUE)),
      c("blue", "white", "red")
    )
    
    column_ha = HeatmapAnnotation(
      Module = rep(module, ncol(top_expr_scaled)),
      col = list(Module = c(setNames(module, module)))
    )
    
    cell_height = 12
    n_rows = nrow(top_expr_scaled)
    n_cols = ncol(top_expr_scaled)
    
    ht = Heatmap(
      top_expr_scaled,
      name = "Expression",
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 10),
      column_names_gp = gpar(fontsize = 8),
      column_title = plot_title,
      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
      top_annotation = column_ha,
      heatmap_legend_param = list(
        title = "Scaled Expression",
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8),
        legend_direction = "horizontal",
        legend_width = unit(8, "cm")
      )
    )
    
    # Return the heatmap object instead of drawing it
    return(ht)
  }
  
  module_heatmaps <- list()
  modules = unique(mergedColors)
  modules = modules[modules != "grey"]
  
  for(module in modules) {
    cat("Processing module:", module, "\n")
    ht = plotModuleTopGenesWithGO(module)
    if(!is.null(ht)) {
      # Store the heatmap object in the list
      module_heatmaps[[module]] <- ht
    }
  }
  
  return(list(
    go_results = all_go_results, 
    summary = summary_df,
    plots = all_plots,
    module_heatmaps = module_heatmaps,
    top_terms = if(exists("top_terms")) top_terms else NULL,
    go_terms = if(exists("summary_df") && nrow(summary_df) > 0) summary_df else NULL,
    module_terms = if(exists("module_terms")) module_terms else NULL
  ))
}

#########################################################
#########################################################
#########################################################

# aggregate on symbols - use this everytwhere so may as well
# be hardcoded in here
geneid2genename <- read.csv("/well/powrie/shared/music_examples/mm39_ensembl107_gene_id2gene_name.tsv",
                            header=TRUE,
                            stringsAsFactors=FALSE,
                            sep="\t")
rownames(geneid2genename) <- geneid2genename$ensemble_gene_id
colnames(geneid2genename) <- c("gene_id", "gene_name")


