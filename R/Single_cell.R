library(dplyr)
library(magrittr) 

generate_boots <- function(celltype, n) {
  dt <- data.frame(cluster = celltype, id = 1:length(celltype))
  index <- do.call(cbind, sapply(1:n, function(x) {
    splits <- dt %>%
      group_by(cluster) %>%
      sample_n(dplyr::n(), replace = TRUE) %>% #sampling with replacement
      ungroup() %>%
      dplyr::select("id")
  }))
  return(index)
}

get_ave_exp <- function(i, myseurat, samples, myident) {
  # 1. Extract the subset of cells
  cell_ids <- samples[, i]
  # 2. Extract and clean metadata
  # Ensure we only take the rows corresponding to our current sample subset
  meta_subset <- myseurat@meta.data[cell_ids, , drop = FALSE]
  # 3. Extract and sanitize the counts matrix
  # We convert to a standard matrix and back to dgCMatrix to strip hidden V5 attributes
  raw_counts <- myseurat@assays$RNA@counts[, cell_ids]
  
  #Clean up the matrix and rownames
  # Force unique names to prevent the LogMap / Duplicate rownames error
  # as.character() ensures we aren't passing any 'Index' objects from Seurat V5
  clean_genes <- make.unique(as.character(rownames(raw_counts)))
  clean_cells <- make.unique(as.character(colnames(raw_counts)))
  #dgCMatrix
  # Re-construct a "clean" sparse matrix
  library(Matrix)
  sample_clean <- Matrix(as.matrix(raw_counts), sparse = TRUE)
  rownames(sample_clean) <- clean_genes
  colnames(sample_clean) <- clean_cells
  # Ensure metadata rownames match the new unique cell names
  rownames(meta_subset) <- clean_cells
  
  # 4. Create the Seurat Object
  SeuratObject <- suppressWarnings(
    CreateSeuratObject(counts = sample_clean, meta.data = meta_subset)
  )
  
  # 5. Process and Calculate Average Expression
  SeuratObject <- NormalizeData(SeuratObject, verbose = FALSE)
  
  # Calculate Average Expression
  # Note: In Seurat V5, AverageExpression returns a list or a standard matrix depending on version
  ave_list <- AverageExpression(SeuratObject, group.by = myident, return.seurat = TRUE)
  # Handle both Seurat V4/V5 return styles for the data slot
  ave <- GetAssayData(ave_list[["RNA"]], layer = "data")
  return(ave)
}

calculate_avg_exp <- function(myseurat,myident,n_bootstrap,seed) {
  set.seed(seed)
  samples=generate_boots(myseurat@meta.data[,myident],n_bootstrap)
  exp <- lapply(1:n_bootstrap,get_ave_exp,myseurat,samples,myident)
  exp <- do.call(cbind, exp)
  return(exp)
}