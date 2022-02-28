#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH -J rf-class_meta
#SBATCH --output=R-%x.%j.o
#SBATCH --error=R-%x.%j.err
#SBATCH --partition=batch-low,batch-medium,batch-high


#======================================================================#
# Use features in table to predict metadata class with random forest model
#
# February 14, 2020			
#======================================================================#

#Set working directory to your R folder 
workdir=/home/see/Wright_Labs/May_2020/Fracking_2019_MT/rf-all
#export R_LIBRARY_PATH=/home/see/R-3.3.3/library
cd $workdir

# Parameters used for 16S rRNA water UOG dataset
dataset=16S_water
input=16S_rRNA_water_ASVs.txt
meta=16S_water_meta-rev-UPDATED-FINAL_rf.txt

mkdir $dataset

source ~/.bashrc
conda activate biom-convert

R --slave --args --source_dir /home/see/R_scripts/R \
-i $input \
-m $meta \
-c HF_Status \
-s SampleID \
-v $dataset-status_stats.txt \
-a $dataset-status_accuracy.csv \
-p $dataset-status_predictors.csv \
-g $dataset-status_gini_plot.pdf \
-q $dataset-status_confusion_matrix.txt \
-t $dataset-status_top_predictors.txt \
-e $dataset-test_d.txt \
-w $dataset-preds2.txt \
-n 20 \
-r 1000 \
< /home/see/R_scripts/rf_classification_cluster_6.10.2021-curves_gen.R



#-i: input table
#-m: metadata file
#-c: class to predict
#-s: name of column with sample ids in table and metadata files, default="SampleID"
#-a: name of csv file containing average accuracy, default="Avg_Percent_Classified_Correct_Value.csv"
#-p: name of csv file containing list of predictors, default="Avg_Mean_Decrease_Gini.csv"
#-g: name of pdf gini predictor importance plot, default="Variable_Importance_plot.pdf"
#-q: name of txt confusion matrix, default="confusion_matrix.txt"
#-t: name of txt for table of best predictors, default="top_predictors.txt"
#-n: number of best predictors to plot in the pdf, default=20
#-r: number of repeats for the random forest model, default=1000
#-v: name of txt for table with model stats