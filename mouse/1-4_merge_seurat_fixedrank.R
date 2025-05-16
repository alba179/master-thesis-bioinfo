# We are going to create a seurat object by reading each matrix from each sample one by one, filtering the empty droplets, and storing the barcodes in new matrices, which are filtered.                                                                                                                                              
# Load libraries
.libPaths("/home/albax/miniforge3/envs/seurat_v4/lib/R/library")

# library(BUSpaRse) # github (devtools)
library(Seurat) #version 4
# library(SeuratWrappers) # github (devtools) 
# library(BSgenome.Mmusculus.UCSC.mm10)
# library(AnnotationHub)
library(zeallot) # For %<-% that unpacks lists in the Python manner   
library(DropletUtils)
library(tidyverse)
library(Matrix)
library(GGally) # For ggpairs
#library(velocyto.R) # github (devtools)
library(SingleR)
library(scales)                   
library(plotly)
library(patchwork)

## set working env and directories to proccess

setwd("./STARalignment/") #directory where we have the matrices (STARAlignment/{SAMPLE}/{RUN}_{SAMPLE}Solo.out/

mydirs <- dir(pattern = "Solo.out",
                  full.names = TRUE,
                  recursive = TRUE,
                  include.dirs = TRUE)

read_matrices <- function(dir=x,
			  spliced_file="Velocyto/raw/spliced.mtx",
			  unspliced_file="Velocyto/raw/unspliced.mtx",
			  gene_file="Gene/raw/matrix.mtx",
			  barcodes="Velocyto/raw/barcodes.txt",
			  genes="Velocyto/raw/genes.txt",
			  gene_barcodes = "Gene/raw/barcodes.txt",
			  gene_genes = "Gene/raw/genes.txt",
			  ...) {

	gene <- ReadMtx(mtx = paste0(dir, "/", gene_file),
                    cells = paste0(dir, "/", gene_barcodes),
                    features = paste0(dir, "/", gene_genes),
                    ...)
	spliced <- ReadMtx(mtx = paste0(dir, "/", spliced_file),
                       cells = paste0(dir, "/", barcodes),
                       features = paste0(dir, "/", genes),
                       ...)
	unspliced <- ReadMtx(mtx = paste0(dir, "/", unspliced_file),
                         cells = paste0(dir, "/", barcodes),
                         features = paste0(dir, "/", genes),
                         ...)

  return(list(gene, spliced, unspliced, dirs = list(gene = paste0(dir, "/", gene_file),
                          spliced = paste0(dir, "/", spliced_file),
                          unspliced = paste0(dir, "/", unspliced_file))))
}

  # Knee plot function
knee_plot <- function(bc_ranks) {
    knee_plt <- tibble(rank = map(bc_ranks, ~ .x[["rank"]]),
                       total = map(bc_ranks, ~ .x[["total"]]),
                       dataset = names(bc_ranks)) %>%
      unnest(cols = c(rank, total)) %>%
      distinct() %>%
      dplyr::filter(total > 0)
    annot <- tibble(inflection = map_dbl(bc_ranks, ~ metadata(.x)[["inflection"]]),
                    rank_cutoff = map_dbl(bc_ranks, ~ max(.x$rank[.x$total >
                                                                  metadata(.x)[["inflection"]]])),
                    dataset = names(bc_ranks))
    p <- ggplot(knee_plt, aes(rank, total, color = dataset)) +
      geom_line() +
      geom_hline(aes(yintercept = inflection, color = dataset),
                 data = annot, linetype = 2) +
      geom_vline(aes(xintercept = rank_cutoff, color = dataset),
                 data = annot, linetype = 2) +
      scale_x_log10() +
      scale_y_log10() +
      labs(x = "Rank", y = "Total UMIs")
    return(p)
  }


# To process each sample (with their matrices)
inflection_points <- list()
knee_plots <- list()
count <- 0
for (x in mydirs) {
  # Read matrices
  count <- count + 1
  cat("\n", count, "Processing sample: ", x, "\n")
  c(gene, spliced, unspliced, dirs) %<-%  read_matrices(dir = x,
						  spliced_file = "Velocyto/raw/spliced.mtx",
						  unspliced_file = "Velocyto/raw/unspliced.mtx",
						  gene_file="Gene/raw/matrix.mtx",
						  barcodes = "Velocyto/raw/barcodes.tsv",
      				genes = "Velocyto/raw/features.tsv",
    					gene_barcodes = "Gene/raw/barcodes.txt",
						  gene_genes = "Gene/raw/genes.txt",
 						  cell.column = 1,
						  feature.column = 2,
						  mtx.transpose = FALSE)
 
  cat("\ndim(gene):", dim(gene),"\n")
  cat("\ndim(spliced):", dim(spliced), "\n")
  cat("\ndim(unspliced)", dim(unspliced), "\n")
  
  cat("\nGene matrix directory: ", dirs$gene, "\n")
  cat("Spliced matrix directory: ", dirs$spliced, "\n")
  cat("Unspliced matrix directory: ", dirs$unspliced, "\n")

  path_names <- strsplit(x, "/")[[1]]
  sample_code <- path_names[2]  # "SIGAD8"
  exp_name <- strsplit(path_names[3], "_")[[1]][1]  # "SLX-17937"
  sample_name <- paste(sample_code, exp_name, sep = "_")

  # Total counts and barcode ranks
  tot_count <- Matrix::colSums(spliced)
  # cat("summary(tot_count):\n", summary(tot_count))

  # to detect emptyDrops
  bc_rank <- barcodeRanks(spliced)
  bc_uns <- barcodeRanks(unspliced)
  bc_gene <- barcodeRanks(gene)

  # Generate knee plot
  knee_plots[[sample_name]]<- knee_plot(list(spliced = bc_rank, unspliced = bc_uns, gene = bc_gene)) +
    coord_flip() +
    ggtitle(sample_name) + 
    theme(plot.title = element_text(hjust = 0.5))

  # Filter barcodes and genes with a FIXED POINT (in our case, 1000). It is recommended to check first the inflection points form the generated plots. 
  cut <- 1000 # change this value as desired, it it sthe UMI count
  bcs_use <- colnames(spliced)[tot_count > cut]
  cat("\ncutoff:",  cut, "\n")
  tot_genes <- Matrix::rowSums(spliced)
  genes_use <- rownames(spliced)[tot_genes > 0]

  sf <- spliced[,bcs_use]
  uf <- unspliced[,bcs_use]
  gn <- gene[,bcs_use]

  cat("\ndim(gene):", dim(gn),"\n")
  cat("\ndim(spliced):", dim(sf), "\n")
  cat("\ndim(unspliced)", dim(uf), "\n")

  # Save filtered matrices
  writeMM(sf, file.path(x, "filtered_spliced.mtx"))
  writeMM(uf, file.path(x, "filtered_unspliced.mtx"))
  writeMM(gn, file.path(x, "filtered_gene.mtx"))

    # Seurat object making (spliced & unspliced):
  print("Making Seurat object...")
  subEsoph <- CreateSeuratObject(gn, assay = "gene")
  print("Merging assays with Seurat object...")
  subEsoph[["spliced"]] <- CreateAssayObject(sf) # change spliced by sf if doing intersect (conservative approach)
  subEsoph[["unspliced"]] <- CreateAssayObject(uf) # change unspliced by uf if doing intersect (conservative approach)

  print("Finished Seurat object for the sample.")


  inflection_points[[sample_name]] <- metadata(bc_rank)$inflection

  colnames_subEsoph <- colnames(subEsoph) # cell barcodes
  new_cellnames <- paste0(colnames_subEsoph, "-", sample_name) # "-" para separar el barcode de nuestro nombre, y hacer los barcodes únicos para cada sample
  subEsoph <- RenameCells(subEsoph, new.names = new_cellnames)
  
  cat("\ndim(subEsoph)", dim(subEsoph), "\n")

  # Merge Seurat objects:
  
  if (count == 1) {
    print("Created Seurat object")
    Esoph <- subEsoph  # make Merged object
  } else {
    print("Merging Seurat objects")
    Esoph <- merge(Esoph, subEsoph)
    cat("\ndim(Esoph)", dim(Esoph), "\n")

  }

  rm(subEsoph)
  gc() # super importante para ir liberando ram
  cat("Finished processing:", x, "\n##################################  \n################################## \n")

}

print(as.list(inflection_points))

# To view the plots combined in a single one: 
combined_plots <- wrap_plots(knee_plots) + plot_layout(ncol = 2)
combined_plots
ggsave("/home/albax/mcGinn_2021/results/combined_plot.png", combined_plots, dpi = 300,  width = 20, height = 45)


# Adding metadata from an exiting csv file. We must create this csv file with the information we want, in each column: 
# sample_name batch_ID  Sequencing_ID (...)
# SIGAD8_SLX-17937  Batch1  SLX-17937
# (...)
cat("Adding metadata to seurat object...")

metadata <- read.csv("~/mcGinn_2021/metadata/metadata_seurat.csv", header = TRUE, sep = ";")

# add Sample metadata to the seurat_obj[[]] slot, ONE by ONE!
# Esoph[['Library_ID']] <- metadata$Library_ID[match(rownames(Esoph@meta.data), metadata$Sample_name)]

process_barcode <- function(seurat_obj, metadata) {
  # Extract Sample_name from cell barcode
  sample <- sapply(strsplit(rownames(seurat_obj@meta.data), "-"), function(x) paste(x[2], x[3], sep = "-"))

  # for each column from the metadata, add it to seurat object with the same name
  for (col in colnames(metadata)) {
    metadata[[col]] <- as.factor(metadata[[col]])
    seurat_obj[[col]] <- NA
    seurat_obj[[col]] <- metadata[[col]][match(sample, metadata$Sample_name)]
  }
  return(seurat_obj)
}

Esoph <- process_barcode(Esoph, metadata)
Esoph

save.image(file = "~/mcGinn_2021/results/seurat_object_filtfixed_mm10.RData")

# Save the RDS file: 
saveRDS(Esoph, file = "~/mcGinn_2021/results/esoph_star_filtfixed_mm10.rds")
cat("Finished creating seurat object for gene, sf and uf.")
cat("Finished adding metadata to seurat object...")

# To parallelize the adding of meta.data 
# suppressMessages(library(future))
# suppressMessages(library(future.apply))
# suppressMessages(library(BiocParallel))
# register(MulticoreParam(workers = 1))  # 2 cpus
# plan(sequential)


# in case you want parallelization: 
# Esoph <- future.apply(process_barcode(Esoph, metadata))
