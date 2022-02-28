
#Lines 5-17 taken from Qiime1 CSS.r script
# Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, Fierer N, Peña AG, Goodrich JK, Gordon JI, Huttley GA, Kelley ST, Knights D, Koenig JE, Ley RE, Lozupone CA, McDonald D, Muegge BD, Pirrung M, Reeder J, Sevinsky JR, Turnbaugh PJ, Walters WA, Widmann J, Yatsunenko T, Zaneveld J, Knight R. 2010. QIIME allows analysis of high-throughput community sequencing data. Nat Methods 7:335-336.

args <- commandArgs(trailingOnly=TRUE)
if(!is.element('--source_dir', args)){
  stop("\n\nPlease use '--source_dir' to specify the R source code directory.\n\n")
}
sourcedir <- args[which(args == '--source_dir') + 1]
source(sprintf('%s/util.r',sourcedir))

args <- commandArgs(trailingOnly=TRUE)
if(!is.element('--source_dir', args)){
  stop("\n\nPlease use '--source_dir' to specify the R source code directory.\n\n")
}
sourcedir <- args[which(args == '--source_dir') + 1]
source(sprintf('%s/util.r',sourcedir))


set_lib_paths("/home/see/miniconda2/envs/R5/lib/R/library") #Function defined in util.r, from https://milesmcbain.xyz/hacking-r-library-paths/



###############################################################################################################
#Set the working directory to the folder with your table and metadata

#setwd("C:/Users/Jeremy/Desktop/picrust/tarm_scripts_3.15/11.7.19/alpha_collate_average/")
###############################################################################################################



# Set library to JCS's path
#set_lib_paths("/home/see/miniconda2/envs/R6/lib/R/library") #Function defined in util.r, from https://milesmcbain.xyz/hacking-r-library-paths/


library(plyr) #call the library of functions for the plyr package
library(dplyr) #call the library of functions for the dplyr package
library(foreach) #call the library of functions for the foreach package
library(doParallel) #call the library of functions for the doParallel package
library(randomForest) #call the library of functions for the randomForest package
library(tibble) #call the library of functions for the tibble package
library(tidyverse) #call the library of functions for the tidyverse package
library(gridExtra) #call the library of functions for the gridExtra package
library(cowplot) #call the library of functions for the cowplot package
library(optparse) #call the library of functions for the optparse package
library(ggplot2) #call the library of functions for the ggplot2 package 
library(PRROC) #call the library of functions for the PRROC package 




# https://cran.r-project.org/web/packages/randomForest/randomForest.pdf


#setwd("C:/Users/Jeremy/Desktop/picrust/tarm_scripts_3.15/11.16.19/datasets2/")


# make option list and parse command line
option_list <- list(
  make_option(c("--source_dir"), type="character",
              help="Path to R source directory [required]."),
  make_option(c("-i", "--input_path"), type="character",
              help="Input normalized table [required]."),
  make_option(c("-m", "--metadata"), type="character",
              help="Input metadata [required]."),
  make_option(c("-v", "--stats"), type="character",
              help="Name of output table containing model stats [required]."),
  make_option(c("-c", "--class"), type="character",
              help="Name of class to predict [required]."),
  make_option(c("-s", "--shared"), type="character",
              help="Name of column with sampleid's in both data and metadata",
              default="SampleID"),
  make_option(c("-a", "--accuracy"), type="character",
              help="Name of csv output file containing average accuracy",
              default="Avg_Percent_Classified_Correct_Value.csv"),
  make_option(c("-p", "--predictors"), type="character",
              help="Name of csv output file containing all predictors",
              default="Avg_Mean_Decrease_Gini.csv"),
  make_option(c("-g", "--gini"), type="character",
              help="Name of pdf output file containing all predictors",
              default="Variable_Importance_plot.pdf"),
  make_option(c("-q", "--matrix"), type="character",
              help="Name of txt output file containing confusion matrix",
              default="confusion_matrix.txt"),
  make_option(c("-t", "--top"), type="character",
              help="Name of txt output file containing top predictors",
              default="top_predictors.txt"),
  make_option(c("-n", "--number"), type="integer",
              help="Specify top however many predictors to include in top predictors output file",
              default=20),
  make_option(c("-r", "--repeats"), type="integer", default=1000,
              help="Number of times to create a random forest model"),
  make_option(c("-z", "--percent"), type="numeric",
              help="Specify amount of samples (as a proportion, e.g. 0.66) to use for training set",
              default=0.66)
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)



###############################################################################################
# Data should be a normalized table
# This script is designed for OTU/ASV tables with taxonomy; however, any table can be used
# In other words, you can use this for genes, pathways, etc. 
# It'll just give you duplicate columns for some of the outputs

data=opts$input_path
meta=opts$metadata
#jobs=opts$jobs
perc=opts$percent
stats=opts$stats


print(perc)
is.numeric(perc)

n=opts$number #Top n predictors to show in the dot plot at the end

shared=opts$shared #column to merge all three samples by
category=opts$class #name of categorical group to predict

#Number of repeats for random forest model
repeatNumber = opts$repeats #number of times to repeat the random forest model generation process, 100 iterations is typically enough for this kind of data


# Output

accuracy=opts$accuracy
predictors=opts$predictors
predictors_plot=opts$gini
confusion_matrix=opts$matrix
gini_table=opts$top


###############################################################################################


asv_table = read.table(data, header=T, check.names = F, sep = "\t", row.names = 1, quote = "") #input the data
metadata = read.table(meta, header=T, check.names = F, sep = "\t") #input the data

if(!"taxonomy" %in% colnames(asv_table))
{
  asv_table$taxonomy=row.names(asv_table);
}


#Make taxonomy table
tax = asv_table
tax2 = data.frame(tax$taxonomy, tax)
tax2 = tax2 %>% rownames_to_column(shared) #tax row names to column
tax3=tax2[,1:2]

colnames(tax3)=c("Predictor_ID","taxonomy")
#tax3 <- data.frame(lapply(tax3, function(x) {gsub(" ", "", x)})) #Delete spaces in this dataframe

#asv_table$SampleID=asv_table$taxonomy
asv_table$taxonomy=NULL

metadata_1=metadata[,c(shared, category)]
metadata_1=as.data.frame(metadata_1)
colnames(metadata_1)=c(shared,"Category")

t_asv2=t(asv_table)


t_asv2=as.data.frame(t_asv2)

t_asv3 = t_asv2 %>% rownames_to_column(shared)

row.names(metadata) = metadata$SampleID
met_num = dplyr::select_if(metadata, is.numeric)
met_num$SampleID = row.names(met_num)
metadata_2 = merge.data.frame(metadata_1, met_num, by=shared)

metadata2=metadata_2[ which(metadata_2$SampleID %in% t_asv3$SampleID), ] 



merge = merge.data.frame(metadata2, t_asv3, by=shared) #merge/full join the metadata with the hData

merge = merge %>% remove_rownames %>% column_to_rownames(var=shared) #Row names to column



merge$Row.names = NULL

merge = na.omit(merge) #delete rows with NA values

merge = droplevels(merge, exclude = if(anyNA(levels(merge$Category))) NULL else NA)                  


confusion=data.frame("anachronox","anachronox","anachronox")
colnames(confusion)=c("test_Predict","Var2","Freq")

###########################################################################


## create dataframe for prediction probabilities and set placeholder values
preds=data.frame(10,10)
## set column names for prediction dataframe
colnames(preds) = c("Class_1","Class_2")

## create dataframe for test data true classes and set placeholder values
test_d=data.frame("anachronox","anachronox")
## set column names
colnames(test_d) = c("Class","SampleID")

###########################################################################

print(merge$Category)
merge$Category = as.factor(merge$Category)

out = foreach(i=1:repeatNumber, .packages=c('randomForest','plyr','dplyr','caret'), .inorder=F, .multicombine=T) %do% { #the loop containing the random forest process, which will be run a repeatNumber of times to produce a list of data frame outputs


    
  index = createDataPartition(merge$Category, p=perc, list=FALSE) #create a random subset of the data (stratified by class)
  train_data = merge[index, ] #create the training set
  test_data = merge[-index, ] #create the test set
  
  mtry_no = round(sqrt(ncol(merge)-1), digits=0) #calculate the number of predictors to be considered at each split
  
  samplesize = rep(round(perc*(table(train_data$Category)))) #calculate the number of samples to use from each class for the bootstrapping
   
  rf = randomForest(train_data[,-which(names(train_data)=="Category")], train_data$Category, sampsize=samplesize, mtry=mtry_no, importance=TRUE) #run random forest to generate a model using the predictors of training set samples  
  test_Predict = predict(rf, test_data, type="class") #predict the class of each test sample using the random forest model produced above
  
  
  ## make predictions with same model and same test data but output probabilities
  test_Predict2 = predict(rf, test_data, type="prob") #predict the class of each test sample using the random forest model produced above
  
  
  ## Save test data in test_d object
  test_data2 = as.data.frame(test_data)
  test_data2$SampleID = row.names(test_data2)
  colnames(test_data2)[1] = "Class"
  test_data3 = test_data2[,c("Class","SampleID")]
  test_d = rbind(test_d, test_data3)
  
  
  ## Save prediction probability output 
  test_Predict3 = as.data.frame(test_Predict2)
  ## Add probability predictions to preds object
  #print(names(test_Predict3))
  colnames(test_Predict3) = c("Class_1", "Class_2")
  preds = rbind(preds, test_Predict3)
  
  
  test_Conting = as.data.frame(table(test_Predict,test_data$Category)) #make a contingency table of predicted vs actual class of the samples in the test set
  
  confusion=rbind(as.matrix(confusion), as.matrix(test_Conting))
  
  #test_Conting$Repeat_number=repeatNumber
  test_match = test_Conting[test_Conting$test_Predict==test_Conting$Var2,] #make a table of only the correct predictions
  test_CorrectClass = (sum(test_match$Freq)/length(test_data$Category)) #calculate the percent of samples classified correctly (a measure of model performance)
  
  imp = importance(rf) #create a table containing the variable importance measures for the predictors
  imp = as.data.frame(imp) #convert the table of importance values to a data frame
  colnames(imp)[ncol(imp)] = as.character(i) #rename the column containing the Mean Decrease in Gini Index vales as the foreach loop iteration number
  imp$Predictor_ID = rownames(imp) #create a new variable "Predictor_ID" that is the rownames of the importance data frame
  imp = imp[-c(1:(ncol(imp)-2))] #remove the unused importance measures
  test_CorrectClassList = list(test_CorrectClass, 'test_CorrectClass') #create a list that includes the previously calculated correct classification value
  imp = rbind(imp, test_CorrectClassList) #append the list from above to the importance info
  
  
  

} #end of foreach loop

# Confusion Matrix stuff

confusion = as.data.frame(confusion)

confusion = filter(confusion, test_Predict !="anachronox") #remove dummy data


confusion = droplevels(confusion, exclude = if(anyNA(levels(confusion$Var2))) NULL else NA) #remove dummy data from level                  
confusion = droplevels(confusion, exclude = if(anyNA(levels(confusion$test_Predict))) NULL else NA)                  

confusion$Freq = as.numeric(as.character(confusion$Freq))


confusion2 = aggregate(. ~ test_Predict + Var2, data = confusion, sum) #combine repetitions to get simplified confusion matrix

actual_percent = ddply(confusion2, "Var2", summarise, sum=sum(Freq)) 

confusion2 = merge.data.frame(confusion2, actual_percent, by="Var2")
confusion2$`Percent (%)`=100*confusion2$Freq/confusion2$sum

###########################################################################

merge.rec = function(.list, ...){ 
  if(length(.list)==1) return(.list[[1]])
  Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
} #create a function that merges all elements in a list


###########################################################################


print("curves step")

## Remove placeholder values
preds2 = preds[!rowSums(preds >1.1),] #Remove values that were greater than 1.1 (placeholder values) since >1.0 is impossible to achieve
test_d = filter(test_d, Class != "anachronox") #Keep all classes besides the initial placeholder class

c1=levels(merge$Category)[1]
c2=levels(merge$Category)[2]

test_d$Class = gsub(c1, 'Class_1', test_d$Class)
test_d$Class = gsub(c2, 'Class_2', test_d$Class)


# Create dataframes with only the class of interest (one dataframe for Class_2 and one for Class_1)
test_d_pos = test_d[(test_d$Class=="Class_2"),] #Create Class_2 dataframe
test_d_minus = test_d[(test_d$Class=="Class_1"),] #Create Class_1 dataframe

test_d_pos$Label = row.names(test_d_pos) #Make column Label with the row names as its values
preds2_pos = subset(preds2, row.names(preds2) %in% test_d_pos$SampleID) #Take only the predictions for samples in test_d_pos dataframe
preds2_pos$Class = "Class_2" #Set value for Column class to Class_2

test_d_minus$Label = row.names(test_d_minus) #Make column Label with the row names as its values
preds2_minus = subset(preds2, row.names(preds2) %in% test_d_minus$SampleID) #Take only the predictions for samples in test_d_minus dataframe
preds2_minus$Class = "Class_1" #Set value for Column class to Class_1

colnames(preds2_minus)[1] = "Class1"
colnames(preds2_minus)[2] = "Class2"
colnames(preds2_pos)[1] = "Class1"
colnames(preds2_pos)[2] = "Class2"

print("Class 1 = Class_1, Class 2 = Class_2")

# Calculate curve values
roc1 = roc.curve(scores.class0=preds2_pos$Class2, weights.class0=preds2_pos$Class2, 
                 scores.class1 = preds2_pos$Class1, weights.class1 = preds2_pos$Class1) #Calculate ROC AUC value for Class_2
print("Class_2 ROC ") #Print enclosed phrase
print(roc1) #Print enclosed object

roc2 = roc.curve(scores.class0=preds2_minus$Class1, weights.class0=preds2_minus$Class1, 
                 scores.class1 = preds2_minus$Class2, weights.class1 = preds2_minus$Class2) #Calculate ROC AUC value for Class_1
print("Class_1 ROC ") #Print enclosed phrase
print(roc2) #Print enclosed object

pr1 = pr.curve(scores.class0=preds2_pos$Class2, weights.class0=preds2_pos$Class2, 
               scores.class1 = preds2_pos$Class1, weights.class1 = preds2_pos$Class1) #Calculate PR AUC value for Class_2
print("Class_2 pr ") #Print enclosed phrase
print(pr1) #Print enclosed object

pr2 = pr.curve(scores.class0=preds2_minus$Class1, weights.class0=preds2_minus$Class1, 
               scores.class1 = preds2_minus$Class2, weights.class1 = preds2_minus$Class2) #Calculate PR AUC value for Class_1
print("Class_1 pr ") #Print enclosed phrase
print(pr2) #Print enclosed object




###########################################################################


#out2=out


out = merge.rec(out, by.x="Predictor_ID", by.y="Predictor_ID", all=T) #merge the output of the foreach loop into a single table by the Predictor_ID column they have in common

test_correctClass = filter(out, Predictor_ID=="test_CorrectClass") #extract the row containing the percent classified correct values
test_correctClassAvg = rowMeans(test_correctClass[-1]) #calculate the average percent classified correctly across all of the iterations of the random forest model
write.csv(test_correctClassAvg, file=accuracy) #write an output file containing the average Percent Classified Correctly value

out = filter(out, Predictor_ID!="test_CorrectClass") #remove the row containing the percent classified correct values from the table
out$GiniAvg = rowMeans(out[-1]) #create a new variable "GiniAvg" that is the mean of all MeanDecreseGini values in each row
out = merge.data.frame(out, tax3, by="Predictor_ID", all=T)
out = out[order(out$GiniAvg, decreasing=TRUE),] #sort the rows by decreasing "GiniAvg" values so the most important predictors are at the top of the table
out = out[,c("Predictor_ID", "GiniAvg", "taxonomy")] #keep only the "Predictor_ID" and "GiniAvg" columns in the table

row.names(out)=c(1:nrow(out)) #Set row names from 1 to top whatever
out$Predictor_Rank=row.names(out) #Send row names to column 'Predictor_Rank'

write.csv(out, file=predictors) #write an output file containing the "GiniAvg" values for all predictors

write.table(confusion2, confusion_matrix, sep = "\t", row.names = F)


#oops = read.csv("Avg_Mean_Decrease_Gini.csv")
## changed out to "oops" below and it to [1:50] to show top 50 predictors

top30 = out[1:n,] #make a new table containing only the top 30 predictors with the highest Mean Decrease in Gini Average
top30 = top30[order(top30$GiniAvg, decreasing=F),] #sort the rows by increasing "GiniAvg" so that the dotchart plots the most important predictors at the top of the graph
pdf(file=predictors_plot) #creates a pdf file that will contain the dotchart below
dotchart(top30$GiniAvg, labels=top30$Predictor_ID, main="Variable Importance", xlab="Mean Decrease in Gini Index", pch=19, xlim=c(0,max(top30$GiniAvg)), cex=0.7) #produce a dotchart of the top 30 most important variables
dev.off()


##################################################################################


## Save PR and ROC AUC values in dataframe
## create dataframe for test data true classes and set placeholder values
aucVal=data.frame("anachronox",5,5,5,5,7,7,7)
## set column names
colnames(aucVal) = c("Dataset","Class_2_PR","Class_2_ROC","Class_2_Accuracy",
                     "Class_1_PR","Class_1_ROC","Class_1_Accuracy",
                     "Overall_Accuracy")

aucVal$Dataset = data

aucVal$`Class_2_PR` = pr1$auc.integral
aucVal$`Class_2_ROC` = roc1$auc
aucVal$`Class_2_Accuracy` = confusion2[4,5]/100

aucVal$`Class_1_PR` = pr2$auc.integral
aucVal$`Class_1_ROC` = roc2$auc
aucVal$`Class_1_Accuracy` = confusion2[1,5]/100

aucVal$Overall_Accuracy = test_correctClassAvg

write.table(aucVal, file=stats, quote=F, row.names = F, sep = "\t")




##################################################################################




# Make nicer dot plot


taxmat=top30

taxmat = separate(taxmat, taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") #It'll fill the levels you don't have with NA

#Delete prefixes (Qiime1) https://stackoverflow.com/questions/29271549/replace-all-occurrences-of-a-string-in-a-data-frame Tim Biegeleisen
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("k__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("p__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("c__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("o__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("f__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("g__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("s__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub(" ", "", x)}))

#Delete prefixes (Qiime2) 
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("D_0__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("D_1__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("D_2__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("D_3__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("D_4__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("D_5__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub("D_6__", "", x)}))
taxmat <- data.frame(lapply(taxmat, function(x) {gsub(" ", "", x)}))





taxmat$GiniAvg=as.numeric(as.character(taxmat$GiniAvg))
taxmat$Predictor_Rank=as.numeric(as.character(taxmat$Predictor_Rank))


#View(taxmat)


#taxmat = transform(taxmat, NodePurity = as.numeric(NodePurity))

is.numeric(taxmat$GiniAvg)
is.numeric(taxmat$Predictor_Rank)


taxmat$plot_order=rev(taxmat$Predictor_Rank) #It plots the smallest numbers at the bottom, so we simply need to reverse the order of the ranks to get it the way we want

write.table(taxmat, gini_table, sep="\t", row.names = F)

print("Done") #print a congratulations message when all lines of code ran

paste("Class_1=", levels(merge$Category)[1])
paste("Class_2=", levels(merge$Category)[2])

