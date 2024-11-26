Population structure
================
Margaux Lefebvre
2024-11-26

# PCA

## LD-pruning

Version: angsd v0.940, python v3.8.12, R v4.2.

NB : I removed the outgroup, since it doesn’t make sense to keep *P.
vivax-like* samples.

``` bash
# Create Beagle input for the pruning
angsd -b list_bams_noPL.txt -ref P.vivax-PvP01_reference_genome.fasta -out PCA_coregenome \
        -rf core_regions.txt \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepthInd 152 -doCounts 1 \
        -GL 2 -doGlf 2 -nThreads 8 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -minMaf 0.05

# Prepare a pos file with the mafs filre
zcat PCA_coregenome.mafs.gz | cut -f 1,2 |  sed 's/:/_/g'| gzip > PCA_coregenome.pos.gz

N_lines=$(zcat PCA_coregenome.pos.gz | wc -l)
let N_pos=N_lines-1 #number of sites

ngsLD \
--geno PCA_coregenome.beagle.gz \
--posH PCA_coregenome.pos.gz \
--probs \
--n_ind 426 \
--n_sites $N_pos \
--max_kb_dist 1 \
--n_threads 8 \
--out PCA_coregenome.ld

python ./scripts/prune_ngsLD.py \
--input PCA_coregenome.ld  \
--max_dist 5000 \
--min_weight 0.5 \
--output PCA_coregenome.snp.unlinked.id
```

Generate an LD-pruned SNP list:

``` r
pruned_position <- as.integer(gsub(paste0("PvP01_[0-1][0-9]_v1:"), "", readLines(paste0("PCA_coregenome.snp.unlinked.id"))))

snp_list <- read.table(paste0("PCA_coregenome.mafs.gz"), stringsAsFactors = F, header = T)[,1:4]

pruned_snp_list <- snp_list[snp_list$position %in% pruned_position, ]
  
write.table(pruned_snp_list, paste0("PCA_coregenome.snp.LDpruned.list"), col.names = F, row.names = F, quote = F, sep = "\t")
```

We keep 247,890 SNPs.

## PCA with PCAngsd

Version: angsd v0.940, PCAngsd v0.98.

``` bash
# Create input
angsd sites index PCA_coregenome.snp.LDpruned.list # mandatory

angsd -b list_bams_noPL.txt -ref P.vivax-PvP01_reference_genome.fasta -out PCA_coregenome_total \
        -rf core_regions.txt \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepthInd 106 -doCounts 1 \
        -GL 2 -doGlf 2 -nThreads 8 -doMajorMinor 1 -doMajorMinor 3 -doMAF 1 -doPost 1 -doIBS 1 -doCov 1 -makeMatrix 1 -sites PCA_coregenome.snp.LDpruned.list
        
# Calculate the PCA with PCAngsd
pcangsd -b PCA_coregenome_total.beagle.gz -o PCA_coregenome.mat
```

# Ancestry plots

Version: PCAngsd v0.98, pong v1.5.

Ancestry plots are inferred with PCAngsd and the input file is the same
as for PCA.

According to [Meisner and Albrechtsen](10.1534/genetics.118.30133), the
best K is determined by 1+ D (the optimal number of principal
components). As presented in Supplementary Fig. 2 in the paper, D would
be equal to 5 (determined by the elbow (broken-stick) method), so K=6.

``` bash
for k in {2..16}
do
 pcangsd -b PCA_coregenome_total.beagle.gz --admix --admix_K $k -o ancestry
done
```

The visualization was done with pong:

``` bash
pong -m file_map.txt -i ind2pop.txt -n country_order.txt -l color.txt 
```

# Tree of the samples (maximum likelihood)

Version: [vcf2phylip v2.0](https://doi.org/10.5281/zenodo.2540861),
python v2.7.5, iqtree v2.0.3, bcftools v1.10.2.

I used *P. vivax-like* as an outgroup.

``` bash
# keep only the core genome
bcftools view VivaxSimium_filtered_final.ploidy2.vcf.gz -R core_genome.txt -o VivaxSimium_core.snpbi_filtered.ploidy2.vcf.gz

# change vcf to phylip file
python ./vcf2phylip/vcf2phylip.py -i VivaxSimium_core.snpbi_filtered.ploidy2.vcf.gz --output-prefix tree_allsamples -o p1537.PL.Cameroon

# ML tree
iqtree -s tree_allsamples.min4.phy -m MFP+ASC -T 24 --prefix tree_allsamples -o p1537.PL.Cameroon -B 1000 -alrt 1000 -st DNA #MFP+ASC = model finder for dataset with only variable sites
```

# IBD Network

## Create the dataset and calculate the IBD

``` bash
mkdir hmm_format
cd hmm_format

bcftools annotate -x INFO,^FORMAT/GT VivaxSimium_filtered_final.ploidy2.vcf.gz | grep -v "##" |cut -d$'\t' -f1-2,10- > temp0

sed 's/0\/0/0/g' temp0 > temp1
sed 's/1\/1/1/g' temp1 > temp2
sed 's/1\/0/-1/g' temp2 > temp3
sed 's/0\/1/-1/g' temp3 > temp4
sed 's/0\/1/-1/g' temp4 > temp5
sed 's/.\/./-1/g' temp5 > temp6
sed 's/PvP01_\([0-9][0-9]*\)_v1/\1/g' temp6 > temp7
sed 's/vPvP01_\([0-9][0-9]*\)_v1/\1/g' temp7 > all_hmm.pf
hmmIBD -i all_hmm.pf -o ../IBD
```

## Read the data and plot the network

``` r
library(readr)
library(igraph)
IBD_df <- read_delim("./IBD.hmm_fract.txt", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)

# Unique nodes
nodes <- sort(unique(c(IBD_df$sample1, IBD_df$sample2)))

# Create an empty symmetric matrix
mat <- matrix(0, nrow = length(nodes), ncol = length(nodes))
rownames(mat) <- colnames(mat) <- nodes

# Fill in the values from the dataset
for (i in 1:nrow(IBD_df)) {
  from <- IBD_df$sample1[i]
  to <- IBD_df$sample2[i]
  value <- IBD_df$fract_sites_IBD[i]
  mat[from, to] <- value
  mat[to, from] <- value
}

# Keep only when IBD > 0.05
mat[mat<0.05] <- 0

# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix( mat, weighted=T, mode="undirected", diag=F)

# Basic chart
plot(network, vertex.size=10, vertex.label=NA) 
E(network)$width <- log(E(network)$weight*100)
plot(network, vertex.size=10, vertex.label=NA) 

# Add colors
list_samples<-as.data.frame(V(network))
list_samples$sample_ID<-rownames(list_samples)
# add meta information
library(tidyverse)
metadata_simium <- read_delim("./metadata_vivaxsimium.tsv", 
    delim = "\t", quote = "\\\"", escape_double = FALSE, 
    trim_ws = TRUE)
metadata_simium<-metadata_simium[,c(2,3,18,17)]
colnames(metadata_simium)<-c("species","host","country","sample_ID")
total_info<-inner_join(list_samples,metadata_simium)

colrs<-c("#d84880","#67d04a","#943fd3","#b4c43f","#635ec7","#cb9f35","#be46ac","#74ba6f","#c53c3b","#62bdae", "#de6a30","#76a5d1","#81552e","#c891cd","#476939","#564f80","#c5ac78","#863a54","#d68a89")

V(network)$color <- colrs[as.factor(total_info$country)]

plot(network, vertex.size=5, vertex.label=NA,edge.curved=0.2, edge.color="#d3d3d3",arrow.mode=0) 

# See which color correspond to which country
legend_col<-unique(cbind(total_info$country,colrs[as.factor(total_info$country)]))
legend_col
```
