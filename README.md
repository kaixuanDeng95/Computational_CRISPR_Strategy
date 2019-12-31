## Crisper
### Requires
- R>=3.5<br>
- R packages: RandomForest,glmnet,e1071 and kernlab<br>
### How to Use
```
Firstly,the Python script kmer.py should be used to get a feature file with a FASTA format input file.
Then the RF model which is stored in Rdata file "RF.model.Rdata" can be used to predict the Z-score of the sequence.
```
### Demonstration
A FASTA file with two DNA sequences is used to predict Z-score.The FASTA file "example.fasta" is tranform to feature file "example_7mer.txt" by python script "kmer.py",then the feature file is used to predict Z-score.
```
>chr6:36634989-36635089
TCTGGCACCCTGCAAGGCCGCATGATGATGCAACAATGCAACAAAAGACAAGCCCGGGCAAGGCCAGCGGGAGCTCTGCCGGCCAGAGTTGCTGATGCGA
>chr6:36635104-36635204
TGGGGAGGGTGTTTCAGGGCTGCAGGGAAGTGGGAGGCCCCAACTGCCCAGGAGGCAAAACTGGCCTCCTGCTCACTCAGCCATGAGCTTTTCTACCCCA
prediction result
> y_pred
4.7719043 0.1287047
```

