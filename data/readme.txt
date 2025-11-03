### PPIs
cd PPIs/
IntAct PSI-MI TAB 2.5 ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip (2021-11-09)
DIP (8 specieses) 2017/02/05 PSI-MI TAB 2.5 https://dip.doe-mbi.ucla.edu/dip/ (2021-11-09)
MINT PSI-MI TAB 2.7 http://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/current/search/query/species:{mouse, human, fruit fly, yeast} (2021-11-09)
BioGRID v4.4.203 PSI-MI TAB 2.5 https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.203/ (2021-11-09)
MatrixDB http://matrixdb.univ-lyon1.fr (2022-03-09)
Skinnider_2021.xlsx # mouse PPIs from Skinnider et al. 2021 (2021-11-08)


### PPIs/PSI-MI/
cd ../../src/PSI-MI/
mi.owl # relationship between MI, https://www.ebi.ac.uk/ols/ontologies/MI (2021-11-10)

### Reactome pathways
cd ../../data/Reactome_pathway
c2.cp.reactome.v7.5.1.symbols.gmt # http://www.gsea-msigdb.org/gsea/downloads.jsp (2022-07-12)
reactome_pathways_from_MsigDB.txt # (src/reactome_pathway.ipynb) (2022-07-13)
UniProt2Reactome_All_Levels.txt # All levels of the pathway hierarchy (2022-12-15) 
reactome_pathway_category.txt # (src/reactome_pathway.ipynb) (2022-12-16)


### Orthology/
## seven_species
species.txt # seven species (2025-04-22)
homologene.data.txt # HID(HomoloGene group id)\tTaxonomy ID\tGene ID\tGene Symbol\tProtein gi\tProtein accession" https://ftp.ncbi.nih.gov/pub/HomoloGene/build68/ (2025-04-22)

treeoflife.species.evolution.tsv # evolution characterization for each species got from Zitnik et al. 2019 (2025-05-09)


### /DNA_sequence
cd ..
seq_download.sh # (2025-06-22)

### /
cd ../data/VTGs
virally-targeted_genes29.txt # PPIs between 29 viruses and human from Tang et al. 2022 (2025-06-20)

### /virus_host
cd ../virus_host/
wget https://www.genome.jp/ftp/db/virushostdb/old/release229/virushostdb.tsv # viral host information downloaded from Virus-Host DB (2025-08-05)
virus_host_eid2, viruses_eid2 # https://eid2.liverpool.ac.uk/OrganismInteractions/ (2025-08-11)
virus_host_summary.txt # ../Analysis/Network_properties/orthologs_with_human/VHI.R (2025-08-14)

### /PoS_vaccine&drug
cd ../PoS_vaccine&drug
PoS_Lo2020.xlsx # Lo et al. 2020. (2025-09-25)
viruses_disease.xlsx # The correspondence between viruses and diseases classified by Lo et al. (2025-09-25)

















