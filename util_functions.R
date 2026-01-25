# Define transformation
transformations <- list(
  none = function(x) x,
  sqrt = function(x) sqrt(x),
  log = function(x) log(x + 1)
)

# === first function ===
# Fisher enrichment function for GO terms in gene clusters
#*@param cluster_genes A vector of gene identifiers in the cluster
#*@param gos A vector of GO term identifiers to test for enrichment
perform_fisher_enrichment <- function(cluster_genes, gos, go_data_long, all_genes) {
  lapply(gos, function(go_term) {
    # Genes annotated with the current GO term
    go_genes <- unique(go_data_long$Gene[go_data_long$GO_Terms == go_term])
    
    a <- length(intersect(cluster_genes, go_genes)) # Genes in cluster AND annotated with GO
    b <- length(setdiff(cluster_genes, go_genes))   # Genes in cluster NOT annotated with GO
    c <- length(setdiff(go_genes, cluster_genes))   # Genes annotated with GO NOT in cluster
    d <- length(setdiff(all_genes, union(cluster_genes, go_genes))) # NOT in cluster and NOT annotated with GO
    
    contingency_table <- matrix(c(a, b, c, d), nrow = 2)
    test <- fisher.test(contingency_table, alternative = "greater",conf.int=T)
    
    data.frame(
      GO_Term = go_term,
      p_value = test$p.value,
      odds_ratio = test$estimate,
      conf_left = test$conf.int[1] #This will always be finite. Though pessimistic, much better than infinite imputation
    )
  }) %>% bind_rows()
}


# === second function ===
#Function to test the enrichment of a given clustering 
#*@param clusters_definitions A data frame with columns 'Gene' and 'Cluster', representing which gene is in which cluster
#*@param go_data_long A data frame with columns 'Gene' and 'GO_Terms', representing gene to GO term annotations
#*@param configuration a descriptive label for the clustering setup (different transformations, correlation methods, cutoff points, or k-values)
test_cluster_enrichment <- function(clusters_definitions, go_data_long, configuration) {
  
  #this is for computing number of clusters with ≥5 genes
  cluster_sizes <- table(clusters_definitions$Cluster)
  cluster_sizes <- cluster_sizes[cluster_sizes >= 5]
  clusters_definitions <- subset(clusters_definitions, Cluster %in% names(cluster_sizes))
  
  # Merge with GO data
  enriched_data <- clusters_definitions %>%
    left_join(go_data_long, by = "Gene") %>%
    filter(!is.na(GO_Terms))
  
  # Filter clusters with ≥5 annotated genes ##this is for enrichment analysis
  go_genes_counts <- enriched_data %>%
    group_by(Cluster) %>%
    summarise(GeneCount = n_distinct(Gene)) %>%
    filter(GeneCount >= 5)
  
  enriched_data <- enriched_data %>%
    filter(Cluster %in% go_genes_counts$Cluster)
  
  # Prepare cluster list
  cluster_list <- split(enriched_data$Gene, enriched_data$Cluster)
  
  # Run enrichment for each cluster
  res <- lapply(names(cluster_list), function(cluster_id) {
    gos <- table(subset(enriched_data,Cluster==cluster_id)$GO_Terms)
    gos <- names(gos[gos>4]) # do enrichment only for GO terms exist in cluster and also appeared more than 4 times, because the rest won't enrich and doing enrichment on them is just useless and time consuming
    # ✅ Call Fisher enrichment function for GO terms in gene clusters
    enrich <- perform_fisher_enrichment(cluster_list[[cluster_id]], gos, go_data_long, all_genes)
    enrich$Cluster <- cluster_id
    enrich
  }) %>% bind_rows()
  
  # Clean and annotate results
  res <- na.omit(res)
  res$adj_pvalue <- p.adjust(res$p_value, method = "BH")
  res$Configuration <- configuration
  
  # Filter significant terms
  filtered_terms <- res %>% filter(adj_pvalue < 0.05)

  # replacing odds ratio with conf.int to solve the infinite problem
  filtered_terms$odds_ratio = filtered_terms$conf_left
  
  # Compute Mean_odds from adjusted filtered_terms (only for significant ones)
  mean_odds <- mean(filtered_terms$odds_ratio, na.rm = TRUE)
  
  # Prepare outputs
  res_summary <- data.frame(
  Configuration = configuration,
  Total_Clusters_With_More_Than_5_Genes = length(unique(clusters_definitions$Cluster)),
  Mean_odds = mean_odds)
  
  # Return both the summary and the full GO results
  list(
  summary = res_summary,
  enrichment = res)
}
