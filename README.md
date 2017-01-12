# findGSE
findGSE is a tool for estimating size of (heterozygous diploid or homozygous) genomes by fitting k-mer frequencies iteratively with a skew normal distribution model, which is written in R.

To use findGSE, one needs to prepare a k-mer histo file generated with short reads, which contains two tab-separated columns. The first column gives frequencies at which k-mers occur in reads, while the second column gives counts of such distinct k-mers. That is, k is also a required information.

Given multiple fastq.gz files, here is a two-step example for counting k-mers with jellyfish (http://www.cbcb.umd.edu/software/jellyfish/):

```R
  zcat *.fastq.gz | jellyfish count /dev/fd/0 -C -o test_21mer -m 21 -t 1 -s 5G
  jellyfish histo -h 3000000 -o test_21mer.histo test_21mer
```

After getting the histo file, supposing findGSE is installed, we can do the following under R environment:

```R
  library("findGSE")
  findGSE(histo="test_21mer.histo", sizek=21, outdir="hom_test_21mer.histo")
```

For more information, two toy examples about GSE for heterozygous and homozygous genomes are provided under findGSE/example/  .
