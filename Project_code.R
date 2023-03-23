##clear the working directory
rm(list=ls())
##libraries
library("ggplot2")
install.packages("ggpubr")
library(ggpubr)
library(readxl)
library(sandwich)

###Reading in data & check for missing values
###Reading in data & check for missing values
###Reading in data & check for missing values

##1.1 read in clinical data 
colorectal<-read.csv(".../ColorectalCancerPatientData.csv")
#check for missing values
summary(colorectal)
# found one missing value among a total of 63 samples 
which(is.na(colorectal),arr.ind = TRUE)
tail(colorectal,5)
#record #63 has missing values in multiple columns
#age, DFS, DFS events, Adj_Radio, Adj_chem are all missing for record #63
#it doesn't make sense to impute the missing values, we will remove the record instead
clinical<-colorectal[-63,]
summary(clinical)

#1.2 Read in genetic data & check for missing values
genetic<-read.csv("/users/yunxing/Desktop/MU-graduate-school/BusinessIntelligence/data/Project/ColorectalCancerGeneExpressionData.csv")
head(genetic)
summary(genetic[c(1:100),])
#look for missing values
which(is.na(genetic),arr.ind = TRUE)
genetic[c(1:5),]
##no missing values are found

### Visuals
### Visuals
### Visuals 

##Histogram of patient age 
hist(colorectal$Age..in.years.,col="orange")

##Pie chart of Duke's stage

#count number of cases each stage
table(colorectal$Dukes.Stage)

df = data.frame(x <- c(16,14,20,12),
                labels <- c(' stage A','stage B','stage C','stage D'))
pie(x,labels)
ggplot(df, aes(x="", y=x, fill=labels)) +geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +theme_void()

## Expression profiles of 8 genes in 5 different patients 
genetic_10<-genetic[c(1:8),]
head(genetic_10)
p1_5<-genetic_10[,c(5,8,13,21,34)]
barplot(as.matrix(p1_5),col=rainbow(8),space=0.8, ylab="expression level of 8 genes")

### Data Management
### Data Management
### Data Management

##. Dimension reduction, check gene expression variability using standard deviation
#create a variable called "std_ev" to add to the file
#this is a row-wise operation 
#need to install the following packages and library
install.packages("matrixStats")
library(matrixStats)
genetic$stdev<-rowSds(as.matrix(genetic[,c(3:64)]))
#make a histogram to visualize gene varability
histogram(genetic$stdev,xlab="Gene expression standard deviation",col="blue",main="Histogram of gene variability")
#select highly variable genes that have a stand deviation greater than 1.0, and put them into the dataframe high_V
#summary of gene variability: high, mid, and low
high_V<-subset(genetic,genetic$stdev>1.0)    
mid_V<-subset(genetic,0.5<=genetic$stdev & genetic$stdev<1)
low_V<-subset(genetic,genetic$stdev<0.5)
summary(high_V)
summary(mid_V)
summary(low_V)


#Dimention reduction, perform correlation analysis on high_V to reveal highly correlated  genes
#since genes labels are rows, patient IDs are columns, the dataframe needs to be transposed first 
library(dplyr)
library(data.table)
high_V_t<-transpose(high_V)
#add column names to the transposed dataframe
colnames(high_V_t)<-high_V$ID_REF
#extract patient IDs from the high_V and add to the front of high_V_t 
ID_REF<-colnames(high_V)
high_V_t<-cbind(ID_REF,high_V_t)
#remove rows 1 and 2, which are old labels from the original data; and last row which are standard deviation 
high_V_n<-high_V_t[-c(1,2,65),]
## now the rows and columns are swapped, however data transpoing have changed the data type to character, we need to convert them back
#convert gene expression columns to numeric,keep the ID_REF column (1st column) as character 
high_V_n<-high_V_n %>% mutate_at(-c(1),as.numeric)
## high_V_n now has 62 rows (patient IDs) and all genetic expression columns are numeric 
##perform correlation analysis among the genes 
 #need to drop the IR_REF column & make a new dataframe high_V_3
#perform correlation analysis on the entire dataframe
high_V_3<-high_V_n[,-1]
summary(high_V_3)
correlation<-cor(high_V_3)
#Find highly correlated genes, with cutoff set at 0.7
#save correlation matrix to a dataframe & select those with correlation coefficients >0.7
cordf <- as.data.frame(as.table(correlation))
corrgene<-subset(cordf, abs(Freq) > 0.7 & abs(Freq) < 1 )
corrgene
# drop redundant genes from high_V_n
drop<-c("213418_at","1553103_a_at","1554007_at","1554741_s","1552349_a_at","1554897_s_at",
        "1553828_at","1554436_a_at","1553995_at","1554899_s_at")
geneselect<-high_V_n[,!(names(high_V_n)%in% drop)]
## The genetic profile data now have 35 columns and 34 genes, the first column is patient ID

#Data managements: combing the two datasets 
#check if the patient ID matches 
gene_PID<-geneselect$ID_REF
clinical_PID<-clinical$ID_REF
setdiff(gene_PID,clinical_PID)
## results showed the patient_IDs match perfectly, now we can move ahead to perform regression analysis
#first combine the two datasets, remove ID_REF from geneselect to remove redundancy
combined_data<-cbind(clinical,geneselect[,-1])
summary(combined_data)
length(combined_data)

### Hypothesis Testing
### Hypothesis Testing
### Hypothesis Testing

## Hypothesis 1:  Dukes stage 
#H0: patient will have similar DFS regardless of Dukes stage
#Ha: patient diagnosed with Dukes stage C & D have lower DFS than Patients at Stage A & B
DFS_early<-subset(combined_data$DFS..in.months.,(combined_data$Dukes.Stage=="A")|(combined_data$Dukes.Stage=="B"))
DFS_adv<-subset(combined_data$DFS..in.months.,(combined_data$Dukes.Stage=="C")|(combined_data$Dukes.Stage=="D")) 
t.test(DFS_early,DFS_adv)
#make a box-plot comparing early stadge with later stage patients
boxplot(DFS_early, DFS_adv,ylab="DFS in months", names=c("Stage A & B", "Stage C & D"),col=c ("blue","green")) 
#Looks like there is an outlier in DFS_adv, remove it and run the testing again
DFS_adv_1<-DFS_adv[-c(7)]
t.test(DFS_early,DFS_adv_1)
boxplot(DFS_early, DFS_adv_1,ylab="DFS in months", names=c("Stage A & B", "Stage C & D"),col=c ("blue","green")) 

## Hypothesis 2: Age 
#H0: patient will have similar DFS regardless of Age
#Ha: patient 60 and older have lower DFS
DFS_all<-combined_data$DFS..in.months.
DFS_60<-subset(combined_data$DFS..in.months.,combined_data$Age..in.years.>=60)
t.test(DFS_all,DFS_60)

## Hypothesis 3 Gender 
#H0: patient will have similar DFS regardless of gender
#Ha: male and female patients have different DFS
DFS_male<-subset(combined_data$DFS..in.months.,combined_data$Gender=="Male")
DFS_female<-subset(combined_data$DFS..in.months.,clinical$Gender=="Female")
t.test(DFS_male,DFS_female)

# Hypothesis 4, Radiation
#H0: Radiation has no effect on survival
#Ha: Radiation makes a difference in survival
DFS_rad<-subset(combined_data$DFS..in.months.,combined_data$Adj_Radio==1)
DFS_no_rad<-subset(combined_data$DFS..in.months.,combined_data$Adj_Radio==0)
t.test(DFS_rad,DFS_no_rad)

# Hypothesis 5, Radiation
#H0: Chemotherapy has no effect on survival
#Ha: Chemotherapy makes a difference in survival
DFS_chem<-subset(combined_data$DFS..in.months.,combined_data$Adj_Chem==1)
DFS_no_chem<-subset(combined_data$DFS..in.months.,combined_data$Adj_Chem==0)
t.test(DFS_chem,DFS_no_chem)


### Regression Analysis
### Regression Analysis 
### Regression Analysis 

##Logistic regression of DFS event as a function of clinical variables
##Logistic regression

##covert categorical variables to numeric 

#convert gender to 1(Male), 0 (Female)
combined_data$Gender<-ifelse(combined_data$Gender=="Male",1,0)

#convert locations to 1 (Right), 2(Rectum), 3(Colon), 4(left)
combined_data$Location<-ifelse(combined_data$Location=="Right",1,
                                 ifelse(combined_data$Location=="Rectum",2,
                                        ifelse(combined_data$Location=="Colon",3,4
                                               )))
  #convert Duke Stages to 1(A),2 (B),3 (C),4 (D)                 
combined_data$Dukes.Stage<-ifelse(combined_data$Dukes.Stage=="A",1,
                            ifelse(combined_data$Dukes.Stage=="B",2,
                            ifelse(combined_data$Dukes.Stage=="C",3,4)))

#models with all predictors
LGR1<-glm(DFS.event~Age..in.years.+Dukes.Stage+Gender+Location+Adj_Radio+
          Adj_Chem, data=combined_data)
#inspect model output
summary(LGR1)
#compute probabilities using predict
pHat<-predict(LGR1,type="response")
#convert probabilities into binary using cutoff=0.5
yHat<-ifelse(pHat>=0.5,1,0)
#proportion of times when yHat matches DFS event
100*mean(combined_data$DFS.event==yHat)

#model with only Duke stage and location
LGR2<-glm(DFS.event~Dukes.Stage+Location, data=combined_data)
#inspect model output
summary(LGR2)
#compute probabilities using predict
pHat<-predict(LGR2,type="response")
#convert probabilities into binary using cutoff=0.5
yHat<-ifelse(pHat>=0.5,1,0)
#proportion of times when yHat matches DFS event
100*mean(combined_data$DFS.event==yHat)

##Linear regression
##Linear regression
#Here I will model DFS in months as a function of genes

#first genes#1-7 in the file

LR1<-lm(DFS..in.months.~`117_at`+`1405_i_at`+`1438_at`+`1552281_at`
        +`1552309_a_at`+`1552348_at`+`1552365_at`, data=combined_data)
summary(LR1)
##note gene 1405_i-at is significant, we will simplify the model using only this gene
LR1_1<-lm(DFS..in.months.~`1405_i_at`, data=combined_data)
summary(LR1_1)
# model using the genes #8-14
LR2<-lm(DFS..in.months.~`1552455_at`+`1552502_s_at`+`1552511_a_at`+`1552519_at`
        +`1552575_a_at`+`1552742_at`+`1552767_a_at`,data=combined_data)
summary(LR2)
##note that gene 1552502_s_at is most significant, we will simplify the model
LR2_1<-lm(DFS..in.months.~`1552502_s_at`,data=combined_data)
summary(LR2_1)

# model using genes #15-21
LR3<-lm(DFS..in.months.~`1552797_s_at`+`1552834_at`+`1552870_s_at`+`1553102_a_at`
        +`1553296_at`+`1553486_a_at`+`1553589_a_at`,data=combined_data)
summary(LR3)
#note gene 1553589_a_at is most significant, we will simplify the model
LR3_1<-lm(DFS..in.months.~`1553589_a_at`,data=combined_data)
summary(LR3_1)

# model using genes #21-27
LR4<-lm(DFS..in.months.~`1553602_at`+`1553613_s_at`+`1553808_a_at`+`1553830_s_at`
        +`1553970_s_at`+`1553994_at`+`1553995_a_at`,data=combined_data)
summary(LR4)
## note that gene 1553602_at is most significant, we will simplfy the model
LR4_1<-lm(DFS..in.months.~`1553602_at`,data=combined_data)
summary(LR4_1)

# model using genes #28-34
LR5<-lm(DFS..in.months.~`1554018_at`+`1554195_a_at`+`1554242_a_at`
        +`1554332_a_at`+`1554474_a_at`+`1554539_a_at`+`1554679_a_at`
        +`1554741_s_at`, data=combined_data)
summary(LR5)
## note that gene 1554018_at and 1554195_a_at are most significant, we will simplfy the model

LR5_1<-lm(DFS..in.months.~`1554018_at`+`1554195_a_at`,data=combined_data)
summary(LR5_1)

#combine all the important predictors
LR_C<-lm(DFS..in.months.~`1554018_at`+`1554195_a_at`+`1553602_at`+`1553589_a_at`
       +`1552502_s_at`+`1405_i_at`,data=combined_data)
summary(LR_C)
#simplify model by dropping variables

LR_C1<-lm(DFS..in.months.~`1554195_a_at`+`1553602_at`+`1553589_a_at`
          +`1552502_s_at`+`1405_i_at`,data=combined_data)
summary(LR_C1)

LR_C2<-lm(DFS..in.months.~`1553602_at`+`1553589_a_at`
          +`1552502_s_at`+`1405_i_at`,data=combined_data)
summary(LR_C2)

LR_C3<-lm(DFS..in.months.~`1553602_at`+`1405_i_at`,data=combined_data)
summary(LR_C3)

#build the model with clinical factors and compare
LR_clinical<-lm(DFS..in.months.~Age..in.years.+Dukes.Stage+Gender+Location+Adj_Radio+
                  Adj_Chem, data=combined_data)
summary(LR_clinical)

## the following models are not included in the report

LR_C4<-lm(DFS..in.months.~`1553602_at`+
          +`1552502_s_at`+`1405_i_at`,data=combined_data)
summary(LR_C4)

# check if these two genes have interactions

LR_C3_inter<-lm(DFS..in.months.~`1553602_at`+`1405_i_at`+`1553602_at`*`1405_i_at`,data=combined_data)
summary(LR_C3_inter)

LR_clinical_1<-lm(DFS..in.months.~Dukes.Stage+Location, data=combined_data)
summary(LR_clinical_1)

#make scatter plot of the two genes with DFS in months, alone with the regression line 

plot(combined_data$`1405_i_at`,combined_data$DFS..in.months.,col="blue",
     ylab="DFS in months", xlab="gene 1405_i_at expression level",main="gene expression & survival")
abline(LR1_1) 
legend("topleft",legend=paste("Adjusted R^2 is", format(summary(LR1_1)$adj.r.squared,digits=3)))

plot(combined_data$`1553602_at`,combined_data$DFS..in.months., col="red",
     ylab="DFS in months",xlab="gene 1553602_at expression level",main="gene expression & survival")
abline(LR4_1)
legend("topright",legend=paste("Adjusted R^2 is", format(summary(LR4_1)$adj.r.squared,digits=3)))





