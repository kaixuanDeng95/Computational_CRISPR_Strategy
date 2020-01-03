## Crisper
### Requires
- R>=3.5<br>
- R packages: RandomForest,glmnet,e1071 and kernlab<br>
### Install
```
git clone git@github.com:kaixuanDeng95/computational-CRISPR-strategy.git
```
### How to Use
```
Firstly,the Python script kmer.py should be used to get a feature file with a FASTA format input file.
Then the trained RF model which is stored in Rdata file "RF.model.Rdata" can be used to predict the Z_scores of the query sequences.
```
### Demonstration
A FASTA file with two DNA sequences is used for demonstrations.The FASTA file "example.fasta" is tranformed to the feature file of "example_7mer.txt" by the python script of "kmer.py". And then the feature file can be put into the trained RF model to obtain their predicted Z_scores.
#### Input File
```
>chr6:36634989-36635089
TCTGGCACCCTGCAAGGCCGCATGATGATGCAACAATGCAACAAAAGACAAGCCCGGGCAAGGCCAGCGGGAGCTCTGCCGGCCAGAGTTGCTGATGCGA
>chr6:36635104-36635204
TGGGGAGGGTGTTTCAGGGCTGCAGGGAAGTGGGAGGCCCCAACTGCCCAGGAGGCAAAACTGGCCTCCTGCTCACTCAGCCATGAGCTTTTCTACCCCA
```
#### Feature File
```
The feature file is a text file with 2 rows and 16384 columns.
```
#### R Script to Predict Z-score
```
library(randomForest)
x=as.matrix(read.table("example_7mer.txt"))
load("RF.model.Rdata")
y_pred=predict(RF.model,x)
```
#### Prediction Result
```
> y_pred
4.7719043 0.1287047
```
