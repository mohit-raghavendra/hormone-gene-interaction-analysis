## SMOTE oversampling and Hormome gene classification

### Dependencies:

imblearn
matplotlib
numpy
tensorflow

### Notebooks for generating samples using SMOTE and running an SVM classifier

Navigate to the hormone-gene-classification/ directory

##### SMOTE_embedding_space.ipynb

Generates source-target gene embeddings using SMOTE.

##### visualise_embeddings.ipynb

Tensorboard tool to visualise and validate the generated embeddings

##### Classifier/hormone_gene_svm_v1.py

SVM classifier to classify positive and negative association of genes with hormones 


## Gene expression network

### Setup:

Download the GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct dataset (the version 7) from the GTEx portal (https://www.gtexportal.org/home/) and save it in root/gene-coexpression-network directory.

### Dependencies:

Python 3.x  
Numpy  
pandas  
networkx  
pickle  
scipy  
statsmodels  

### Notebooks for obtaining the cardinality vertex cover of the gene expression network.

Navigate to the gene-coexpression-network/ directory

##### conflict_graph_vertex_cover.ipynb

The notebook generates the co-expression network, builds a conflict graph out of the co-expression network and finds the vertex cover of the graph.

##### growth_of_emperical_FDR.ipynb

The notebook analyses the variation of emperical FDRs against increasing the alpha (cutoff) for p values after running multiple testing.
