# Master's thesis in **Bioinformatics and Computational Biology**

This is the repository for my Master's Thesis in **Bioinformatics and Computational Biology** at the Universidad Autónoma de Madrid (UAM), course 2024/25. It contains the code used. 

**👩‍💻 Author**: **Alba Méndez Alejandre**  

## 🧬 Brief description
This work applies variant calling from **scRNA-seq** to link genetic mutations to cellular phenotypes, using a **customized pipeline** on mouse and human esophageal data. It highlights the importance of experimental design and filtering for reliable mutation detection.

## 🎯 Objectives
The main objective of this work is **to assess whether single-cell transcriptomics is suitable for reliable somatic variant detection in mouse esophageal epithelium, and to associate genotype-to-phenotype relations using single cell RNA-seq.**

- Apply and customize variant caller for scRNA-seq in mouse and human esophagus.
- Characterize the diversity of transcriptomic states in the mouse esophageal epithelium. 
- Map variants onto UMAP embeddings. 

## 🛠️ Tools & Technologies used
- (variant caller)
- `Seurat`, `scVelo`, `Slingshot` (transcriptomic analysis & trajectory inference)
- `VEP` (variant annotation)

## 📊 Data
We obtained the HCA data from previous work (González-Menéndez, 2024), originally taken from: 
> Madissoon, E., Wilbrey-Clark, A., Miragaia, R.J. et al. scRNA-seq assessment of the human lung, spleen, and esophagus tissue stability after cold preservation. Genome Biol 21, 1 (2020). https://doi.org/10.1186/s13059-019-1906-x

The mouse data was obtained from a publicly available dataset.

## 📁 Repository structure
The repository holds independent scripts for each dataset: 
```
├── docs/ # quarto book
│
├── Introduction.qmd # installations, brief explanation of configurations used
│
├── human/ # 🧍Human analysis
│ ├── 1_Inspection.qmd # preliminary inspection of the dataset
│ ├── 2_GeneExpression.qmd # calculate average gene expression for sets of genes
│ ├── 3_VennDiagrams.qmd # obtain venn diagrams for sets of genes or mutated genes of interest
│ ├── 4_UMAP_mapping.qmd # map mutated cells in the umap
│ └── mut_clones_analysis_hca.qmd # modifying seurat_obj@meta.data to add clones
│
├── mouse/ # 🐭 Mouse analysis
│ ├── fastQC/
│ │   ├── run1/ # multiqc report for the first sequencing run
│ │   └── run2/ # multiqc report for the second sequencing run
│ ├── 1-4_merge_seurat_fixedrank.R # script to filter out droplets and doublets from the matrices
│ ├── 1_DataProcessing.qmd # notebook which contains the bash scripts used to download the files, and oerform the alignment with STARsolo.
│ ├── 2_ClusteringCellAnnotation.qmd # seurat pipeline to cluster the cells and annotate them
│ ├── 3_Velocity_inference.qmd # infer velocity and pseudotime trajectory
│ ├── 4_VariantCalling.qmd  # steps to perform variant calling and pre-filtering processing, as well as TSV to VCF conversion
│ ├── 5_AnnotationVariants.qmd # exploration of the VEP output
│ └── 6_GOAnalysis.qmd # GO analysis of the mutated genes
│
└──  scripts/ # other general scripts
  └──  annotate_vep.sh # automated annotation with vep

```