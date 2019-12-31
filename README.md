# Crisper
### Requires
- R>=3.5<br>
- R packages:RandomForest,glmnet,e1071 and kernlab<br>
### How to Use
```
Firstly,the Python script kmer.py should be used to get a feature file with a FASTA format input file.Then the Rdata file RF.model.Rdata
can be used to predict the Z-score of the sequence.<br>
```
### Demonstration

The input file must be DNA sequences with a length of 100bp and a format in FASTA.
>chr6:36634989-36635089
TCTGGCACCCTGCAAGGCCGCATGATGATGCAACAATGCAACAAAAGACAAGCCCGGGCAAGGCCAGCGGGAGCTCTGCCGGCCAGAGTTGCTGATGCGA
>chr6:36635104-36635204
TGGGGAGGGTGTTTCAGGGCTGCAGGGAAGTGGGAGGCCCCAACTGCCCAGGAGGCAAAACTGGCCTCCTGCTCACTCAGCCATGAGCTTTTCTACCCCA
