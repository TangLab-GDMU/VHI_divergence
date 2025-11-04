#!/usr/bin/env bash

cd ..
data=$PWD/data
intermediate=$PWD/intermediate
result=$PWD/result
src=$PWD/src

### CDS FASTA files download 
mkdir -p $data/CDS_sequence
cd $data

# genome (2025-06-20)
#wget https://ftp.ensembl.org/pub/release-114/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz
#wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
#wget https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
#wget https://ftp.ensembl.org/pub/release-114/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz
#wget https://ftp.ensembl.org/pub/release-114/fasta/rattus_norvegicus/dna/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz
#wget https://ftp.ensembl.org/pub/release-114/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz

# annotation (2025-06-20) 
#wget https://ftp.ensembl.org/pub/release-114/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.114.gtf.gz
#wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
#wget https://ftp.ensembl.org/pub/release-114/gtf/mus_musculus/Mus_musculus.GRCm39.114.gtf.gz
#wget https://ftp.ensembl.org/pub/release-114/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.114.gtf.gz
#wget https://ftp.ensembl.org/pub/release-114/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.114.gtf.gz
#wget https://ftp.ensembl.org/pub/release-114/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.114.gtf.gz

# CDS
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
wget https://ftp.ensembl.org/pub/release-114/fasta/mus_musculus/cds/Mus_musculus.GRCm39.cds.all.fa.gz
wget https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz
wget https://ftp.ensembl.org/pub/release-114/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.54.cds.all.fa.gz
wget https://ftp.ensembl.org/pub/release-114/fasta/rattus_norvegicus/cds/Rattus_norvegicus.GRCr8.cds.all.fa.gz
wget https://ftp.ensembl.org/pub/release-114/fasta/caenorhabditis_elegans/cds/Caenorhabditis_elegans.WBcel235.cds.all.fa.gz
wget https://download.xenbase.org/xenbase/Genomics/Sequences/xlaevisMRNA.fasta.gz

### Extraction of the longest transcript for each gene (2025-06-21)

seqkit seq -w 0 Homo_sapiens.GRCh38.cds.all.fa | awk '/^>/ {if (seq) print gene "\t" seq; gene = $0; sub(/.*gene_symbol:/, "", gene); sub(/ .*/, "", gene); seq = ""; next} {seq = seq $0} END {print gene "\t" seq}' | sort -u | awk -F '\t' -v ORF='\t' '{print ">" $0, length($2)}' | sort -k1,1 -k3,3nr | awk '!seen[$1]++' | awk -F ' ' '{print $1,$2}' | tr ' ' '\n' > CDS_sequence/Hs.cds.fa

seqkit seq -w 0 Caenorhabditis_elegans.WBcel235.cds.all.fa | awk '/^>/ {if (seq) print gene "\t" seq; gene = $0; sub(/.*gene_symbol:/, "", gene); sub(/ .*/, "", gene); seq = ""; next} {seq = seq $0} END {print gene "\t" seq}' | sort -u | awk -F '\t' -v ORF='\t' '{print ">" $0, length($2)}' | sort -k1,1 -k3,3nr | awk '!seen[$1]++' | awk -F ' ' '{print $1,$2}' | tr ' ' '\n' > CDS_sequence/Ce.cds.fa

seqkit seq -w 0 Drosophila_melanogaster.BDGP6.54.cds.all.fa | awk '/^>/ {if (seq) print gene "\t" seq; gene = $0; sub(/.*gene_symbol:/, "", gene); sub(/ .*/, "", gene); seq = ""; next} {seq = seq $0} END {print gene "\t" seq}' | sort -u | awk -F '\t' -v ORF='\t' '{print ">" $0, length($2)}' | sort -k1,1 -k3,3nr | awk '!seen[$1]++' | awk -F ' ' '{print $1,$2}' | tr ' ' '\n' > CDS_sequence/Dm.cds.fa

seqkit seq -w 0 Mus_musculus.GRCm39.cds.all.fa | awk '/^>/ {if (seq) print gene "\t" seq; gene = $0; sub(/.*gene_symbol:/, "", gene); sub(/ .*/, "", gene); seq = ""; next} {seq = seq $0} END {print gene "\t" seq}' | sort -u | awk -F '\t' -v ORF='\t' '{print ">" $0, length($2)}' | sort -k1,1 -k3,3nr | awk '!seen[$1]++' | awk -F ' ' '{print $1,$2}' | tr ' ' '\n' > CDS_sequence/Mm.cds.fa

seqkit seq -w 0 Rattus_norvegicus.GRCr8.cds.all.fa | awk '/^>/ {if (seq) print gene "\t" seq; gene = $0; sub(/.*gene_symbol:/, "", gene); sub(/ .*/, "", gene); seq = ""; next} {seq = seq $0} END {print gene "\t" seq}' | sort -u | awk -F '\t' -v ORF='\t' '{print ">" $0, length($2)}' | sort -k1,1 -k3,3nr | awk '!seen[$1]++' | awk -F ' ' '{print $1,$2}' | tr ' ' '\n' > CDS_sequence/Rn.cds.fa

seqkit seq -w 0 Saccharomyces_cerevisiae.R64-1-1.cds.all.fa | awk '/^>/ {if (seq) print gene "\t" seq; gene = $0; sub(/.*gene_symbol:/, "", gene); sub(/ .*/, "", gene); seq = ""; next} {seq = seq $0} END {print gene "\t" seq}' | sort -u | awk -F '\t' -v ORF='\t' '{print ">" $0, length($2)}' | sort -k1,1 -k3,3nr | awk '!seen[$1]++' | awk -F ' ' '{print $1,$2}' | tr ' ' '\n' > CDS_sequence/Sc.cds.fa

cat xlaevisMRNA.fasta | grep -A1 "GenePageIDs" | awk '/^>/ {if (seq) print gene "\t" seq; gene = $0; sub(/.*GenePageIDs:/, "", gene); sub(/ .*/, "", gene); seq = ""; next} {seq = seq $0}' |  awk -F '\t' -v ORF='\t' '{print ">" $0, length($2)}' | sed -E 's/^(>[^.|]+)(\.[^|]+)?\|[0-9]+/\1/' | sort -k1,1 -k3,3nr | awk '!seen[$1]++' | awk -F ' ' '{print $1,$2}' | tr ' ' '\n' | awk '/^>/ {print} !/^>/ {print toupper($0)}' > CDS_sequence/Xl.cds.fa
