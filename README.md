# findGSE
findGSE is a tool for estimating size of (heterozygous diploid or homozygous) genomes by fitting k-mer frequencies iteratively with a skew normal distribution model, which is written in R ([code](https://github.com/schneebergerlab/findGSE/blob/master/R/findGSE_v1.94.R)).

To use findGSE, one needs to input a k value and a corresponding k-mer histo file generated with short reads, which contains two tab-separated columns. The first column gives frequencies at which k-mers occur in reads, while the second column gives counts of such distinct k-mers ([example](https://github.com/schneebergerlab/findGSE/blob/master/example/homozygous/test_21mer.histo)).

Given multiple fastq.gz files, here is a two-step example for counting k-mers with [jellyfish](http://www.cbcb.umd.edu/software/jellyfish/):

```R
  zcat *.fastq.gz | jellyfish count /dev/fd/0 -C -o test_21mer -m 21 -t 1 -s 5G
  jellyfish histo -h 3000000 -o test_21mer.histo test_21mer
```

After getting the .histo file, supposing findGSE has been installed ([INSTALL](https://github.com/schneebergerlab/findGSE/blob/master/INSTALL)), we can do the following for GSE under R environment:

```R
  library("findGSE")
  findGSE(histo="test_21mer.histo", sizek=21, outdir="hom_test_21mer")
```

Results will be printed like "Genome size estimate for test_21mer.histo: 1498918 bp." 
For more information about estimation, one can check the .txt and .pdf files in the output dir.

Two detailed [toy examples about GSE for heterozygous and homozygous genomes] (https://github.com/schneebergerlab/findGSE/tree/master/example) are provided for playing around.
