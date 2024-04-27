# OOP Project with Scanpy package

### Goal
Our project involves using the Scanpy package and enhancing the 'plotting' aspect to make the graphs dynamic and interactive for the user.
___

### Package used
Scanpy is an open-source Python library specifically designed for the analysis of single-cell sequencing data (scRNA-seq). It provides a range of functionalities for data preprocessing, visualization, analysis, and interpretation of results.
___
### Setup
#### Installation
To install this version of Scanpy package, there are several steps :

1) Clone the repository to the local machine  
         ```[git clone https://github.com/Ju960/OOP_Project_Julia_GOUNIN_Reshma_VASANTE.git](https://github.com/Ju960/OOP_Project_Julia_GOUNIN_Reshma_VASANTE.git)```

2) Once the clone is complete, you will find a folder named <RepositoryName> on you computer. Navigate from to terminal access folder  
         ``` cd OOP_Project_Julia_GOUNIN_Reshma_VASANTE```

3) Read the various read_me files for each enhancement to better understand how they work

#### Usage
After the installation, you can use this version of scanpy package in Python environment. Before running the programs, remember to install the various packages used, such as scanpy.  
         ```pip install  ```

In the folder, we will find 3 files :  
1. a file that uses widgets
2. a file that uses Plotly
3. a file that uses Bokeh

Each .py file contains the various functions for the enhancements and the notebook will allow you to test your files.
The nottebook must be on the same folder as the 3 files.
#### Data used
In the notebook, you will find all the necessary data to construct the matrix plots :  
- Dataset pmbc
- List of genes of interest
- Clusters
- Marker genes

This dataset (gene_list and marker_gene) can be modified according to your needs and research requirements.

This our version of our dataset : 

```#Dataset
pbmc = sc.datasets.pbmc68k_reduced()

#List of genes of interest
gene_list = ['CD79A', 'MS4A1', 'FCER1A', 'CST3', 'FCGR3A', 'GNLY', 'NKG7', 'IGLL1', 'IGJ', 'CD3D']

#Calculation of clusters
sc.tl.leiden(pbmc, key_added='clusters', resolution=0.5)

#Mapping clusters to cell annotations
cluster2annotation = {
    '0': 'Monocytes',
     '1': 'Dendritic',
     '2': 'T-cell',
     '3': 'NK',
     '4': 'B-cell',
     '5': 'Dendritic',
     '6': 'Plasma',
     '7': 'Other',
     '8': 'Dendritic',
}
pbmc.obs['cell type'] = pbmc.obs['clusters'].map(cluster2annotation).astype('category')

#Marker gene set
marker_genes_dict = {
    'B-cell': ['CD79A', 'MS4A1'],
    'Dendritic': ['FCER1A', 'CST3'],
    'Monocytes': ['FCGR3A'],
    'NK': ['GNLY', 'NKG7'],
    'Other': ['IGLL1'],
    'Plasma': ['IGJ'],
    'T-cell': ['CD3D'],
}
```
___
### Authors
**Students of M1 GENIOMHE**  
Julia GOUNIN  
Reshma VASANTE


