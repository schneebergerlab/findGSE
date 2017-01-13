# You can use the following steps 0 to 2 to simulate any histo files, or you can directly go to step 4 for testing findGSE.

# read simulator to install: ftp://ftp.genomics.org.cn/pub/pIRS/ReadMe.txt

cd example/homozygous/

# step 0. simulate a genome: 3 sequences of 500000 bp long; A/C/G/T equal at 25 percent. In total 1.5 Mb  saved in file random_sequences_set.fa.

perl ../random_dna_strings.pl

# step 1. simulate reads

pirs simulate -i random_sequences_set.fa -m 250 -l 99 -x 50 -v 10 -e 0.01 > sim.log

# step 2. count k-mers with jellyfish -- caution: require about 3.5Gb RAM

zcat *.fq.gz | jellyfish count /dev/fd/0 -C -o test_21mer -m 21 -t 1 -s 1G
jellyfish histo -h 3000000 -o test_21mer.histo test_21mer

# step 3. estimate homozygous genome size (under R)
R
library("findGSE")
findGSE(histo="test_21mer.histo", sizek=21, outdir="hom_test_21mer")
