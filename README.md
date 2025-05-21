# Master's thesis in **Bioinformatics and Computational Biology**

This is the repository for my Master's Thesis in **Bioinformatics and Computational Biology** at the Universidad Autónoma de Madrid (UAM), course 2024/25. It contains the code used. 

**👩‍💻 Author**: **Alba Méndez Alejandre**  

## 🧬 Brief description
This thesis applies variant calling from **scRNA-seq** to link genetic mutations to cellular phenotypes, using a **customized pipeline** on mouse and human esophageal data. It highlights the importance of experimental design and filtering for reliable mutation detection.

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

Mouse data was taken from:
> McGinn, J., Hallou, A., Han, S., Krizic, K., Ulyanchenko, S., Iglesias-Bartolome, R., England, F. J., Verstreken, C., Chalut, K. J., Jensen, K. B., Simons, B. D., & Alcolea, M. P. (2021). A biomechanical switch regulates the transition towards homeostasis in oesophageal epithelium. Nature Cell Biology, 23(5), 511–525. https://doi.org/10.1038/s41556-021-00679-

## 📁 Repository structure
The repository holds independent scripts for each dataset: 
```
└── human/ # 🧍Human analysis
├── 
├── 
└── 
│
├── mouse/ # 🐭 Mouse analysis
│ ├── 
│ ├── 
│ └── 
│
├── scripts/ # other general scripts
│ └──  vep_annotation.sh # automated annotation with vep

```



