# TAPVC Rscript 
 Our script simultaneously trains and validates the TAPVC data by four machine learning methods (LR/RF/GBM/SVM).
# Machining lerning
 Rscript machine.learing.R $oudir train.featuresMerge.exp train.sam.list pred.featuresMerge.exp pred.sam.list case control 
# train2pred.R
 A custom function that uses the training model for prediction.
# Train.featuresMerge.exp
 The first column is the gene name, the other columns are the sample names of the training set, and the expression is normalized by TPM.
# Pred.featuresMerge.exp
 The first column is the gene name, the other columns are the sample names of the validation set, and the expression is normalized by TPM.
# Train.sam.group
 The group information corresponding to the sample name of the training set (0 represents control/1 represents case).
# Pred.sam.group
 The group information corresponding to the sample name of the validation set (0 represents control/1 represents case).
