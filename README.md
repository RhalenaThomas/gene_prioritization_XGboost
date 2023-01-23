# Gene prioritization of Parkinson's disease GWAS loci

Collection of script used to extract multi-modality data from eQTL, single-cell RNAseq data, create neighborhood score, extract distance-based information, run machine learning algorithm. 

# 1. Create the training and testing set.
Run the scripts in the XGBoost folder to create the datasets.

# 2. Generate the machine learning model
Run the Jupyter Notebook- xgboost_v10_imbalanced_v2.ipynb in the XGBoost folder to create the machine learning model.

# 3. Extract the important features
Run the Jupyter Notebook- extract_model_info.ipynb in the XGBoost folder to the important feature from the machine learning model.
