# Gene prioritization of Parkinson's disease GWAS loci

Collection of script used to extract multi-modality data from eQTL, single-cell RNAseq data, create neighborhood score, extract distance-based information, run machine learning algorithm. 

This model is used in the manuscript:
Eric Yu, Roxanne Larivière, Rhalena A Thomas, Lang Liu, Konstantin Senkevich, Shady Rahayel, Jean-François Trempe, Edward A Fon, Ziv Gan-Or, Machine learning nominates the inositol pathway and novel genes in Parkinson’s disease, Brain, Volume 147, Issue 3, March 2024, Pages 887–899, https://doi.org/10.1093/brain/awad345


# Usage

# 1. Create the training and testing set.
1. Collect all the input data from various databases to collect SNP and gene features
2. Run the R scripts within each datafolder to generate the input txt files.
3. Run the R scripts in the XGBoost folder to create the datasets formatted for the XGboost model.

# 2. Generate the machine learning model
Run the Jupyter Notebook- xgboost_v9_add_cv_imbalanced_v2.ipynb in the XGBoost folder to perform feature selection.
Run the Jupyter Notebook- xgboost_v10_imbalanced_v2.ipynb in the XGBoost folder to create the machine learning model.

# 3. Extract the important features
Run the Jupyter Notebook- extract_model_info.ipynb in the XGBoost folder to the important feature from the machine learning model.
