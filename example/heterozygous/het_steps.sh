# You can use the following steps 0.1 to 3 to simulate any histo files, or you can directly go to step 4 for testing findGSE.

# read simulator to install: ftp://ftp.genomics.org.cn/pub/pIRS/ReadMe.txt

# step 0.1. if not done, compile mutate_dna ready for mutating DNA sequence

cd example/mutate_dna/
make

# step 0.2. get into the folder

cd ../heterozygous/

# step 0.3. simulate a genome: 3 sequences of 500000 bp long; A/C/G/T equal at 25 percent. In total 1.5 Mb.

perl ../random_dna_strings.pl

# step 1. mutate genome

../mutate_dna/mutate_dna random_sequences_set.fa 0.005 mutated_

# step 2. simulate reads

snprate=0.005
cov=80
err=0.001
pirs simulate -i random_sequences_set.fa -I mutated_.rate${snprate}.fasta -m 250 -l 99 -x ${cov} -v 10 -e ${err} -o snprate${snprate}_cov${cov} > snprate${snprate}_cov${cov}.log

# step 3. count k-mers with jellyfish -- caution: require about 3.5Gb RAM

zcat *.fq.gz | jellyfish count /dev/fd/0 -C -o test_21mer -m 21 -t 1 -s 1G
jellyfish histo -h 3000000 -o test_21mer.histo test_21mer

# step 4. estimate homozygous genome size
R
library("findGSE")
findGSE(histo="test_21mer.histo", sizek=21, outdir="het_test_21mer", exp_hom=80)
