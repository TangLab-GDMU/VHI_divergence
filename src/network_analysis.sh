#!/usr/bin/env bash

cd ..
data=$PWD/data
intermediate=$PWD/intermediate
result=$PWD/result
src=$PWD/src



#######################
### sequence similarity
### sequence download
./sequence_download.sh

### pairwise sequences extraction
mkdir -p $intermediate/sequence_similarity
python $src/pairwise_sequences_extraction.py \ 
    -oh $data/orthology/seven_species/ortho1to1_with_human.txt \
    -ifa $data/DNA_sequence/CDS_sequence/ \
    -ofa $intermediate/sequence_similarity/ 

### pairwise sequences alignment, and alignment trimming
mkdir -p $intermediate/sequence_similarity/trimmed_aligned_fasta
cd $intermediate/sequence_similarity/
for f in *.fasta
do 
    echo $f
    trimal -in $f -out trimmed_aligned_fasta/trimmed_$f -automated1
    rm $f
done

### pairwise sequences similarity (1-distance)
cd $src/
Rscript --vanilla $src/sequence_similarity.R

tar -czf $intermediate/trimmed_aligned_fasta.tar.gz $intermediate/trimmed_aligned_fasta/
rm -rf $intermediate/trimmed_aligned_fasta/


####################
### degree exponents
cd $result/integrated_ppi

for d in d in *
do
    echo $d
    tail -n +2 $d/integrated_ppi.txt | cut -f1,2 > $d/integrated_ppi2.txt
done

for d in d in *
do
    echo $d
    python $src/degree_exponent_get.py -f $d/integrated_ppi2.txt
done > degree_exponent.txt


##################################
### PPIs sampling (random network)
mkdir -p $intermediate/random_network

python $src/random_network_get.py \
    -d $result/ \
    -s 'H. sapiens' \
    -se 2914 \
    -o $intermediate/random_network |
    tee -a $result/random_network_size.txt


########################################
### interaction conservation coefficient

python ICC.py \
    -d $result/ \
    -s 'H. sapiens' \
    -g $data/orthology/gene_groups.txt \
    -oh $data/orthology/ortho1to1_with_human.txt \
    -o $intermediate/random_network

mkdir -p $result/ICC/
python ICC2.py \
    -i $intermediate/random_network \
    -s 'H. sapiens' \
    -o $result/ICC/


#######################
### LCC identification

### for original interactomes of seven species
num_cores=29

for species in Ce  Dm  Hs Mm  Rn  Sc  Xl
do
    echo $species
    mkdir -p $intermediate/viral_module/$species/intermediate/
    
    parallel -u -j $num_cores --bar \
        python viral_module.py \
            -v $data/VTGs/$species/VTGs/{}.txt \
            -i $result/integrated_ppi_$species.txt \
            -o1 $intermediate/viral_module/$species/{}_lcc_statistic.txt \
            -o2 $intermediate/viral_module/$species/{}_expected_lcc_size.txt \
            ::: `ls $data/VTGs/$species | sed 's/.txt//'`
done

### for 1000 randomlized interactomes of human
for species in Ce  Dm  Mm  Rn  Sc  Xl
do
    echo $species
    
    for i in {0..999}
        do 
            echo edge_list_$i.txt
            sed -i 's/ /\t/' intermediate/random_network/$species/edge_list_$i.txt
            parallel -u -j $num_cores --bar \
                python viral_module.py \
                    -v $data/VTGs/Hs/VTGs/{}.txt \
                    -i intermediate/random_network/$species/edge_list_$i.txt \
                    -o1 $intermediate/viral_module/$species/intermediate/edge_list_"$i"_{}_lcc_statistic.txt \
                    -o2 $intermediate/viral_module/$species/intermediate/edge_list_"$i"_{}_expected_lcc_size.txt \
                    ::: `ls $data/VTGs/Hs/VTGs/ | sed 's/.txt//'`
        done
done


############################
### isolated components (IC)
mkdir -p $result/IC

python IC.py \
    -d $result/ \
    -o $intermediate/random_network/

python IC_orthologs.py \
    -oh $data/orthology/ortho1to1_with_human.txt \
    -o $intermediate/random_network/ \
    -of $result/IC/


#############################
### results and visualization
mkdir -p $result/ICC/GSEA
mkdir -p $result/sequence_similarity/GSEA
mkdir -p $result/viral_module

### 1:1 orthologs (Fig S3)
Rscript --vanilla $src/orthology/homologene.R

### the relationship between ICC and sequence similarity (Fig 1, Fig 2, Fig S4, Fig S5, Fig S6)
Rscript --vanilla $src/ICC_SS.R

### viral modules (Fig 3, Fig S7, Fig S8, Fig S9, Fig S10)
Rscript --vanilla $src/VHI.R

### viral perturation (Fig 4, Fig S11, Fig S12, Fig S13, Fig S14, Fig S15)
Rscript --vanilla $src/network_resilience.R

### rate of phase 1 to approval (Fig 5)
Rscript --vanilla $src/therap_approval.R























