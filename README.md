# Genomic Insights into the Host shifts and Genetic Exchange Between *Plasmodium vivax* and *Plasmodium simium*

This repository is for this paper:

**Genomic insights into host shifts between *Plasmodium vivax* and *Plasmodium simium* in Latin America**,
Margaux J. M. Lefebvre, Fanny Degrugillier, Céline Arnathau, Camila González, Silvia Rondón, Andrés Link, Andrea Chaves, Julio A. Benavides, Aline Alves Scarpellini Campos, Edmilson dos Santos, Rosana Huff, Cláudia Maria Dornelles Silva, Ezequiel Vanderhoeven, Benoit De Thoisy, Michael C. Fontaine, Franck Prugnolle, Virginie Rougeron. *Preprint*. doi: [https://doi.org/10.1101/2024.12.19.629455](https://doi.org/10.1101/2024.12.19.629455)

The languages used are mainly bash and R. For each part, software and version are specified.

# Summary

## [About the dataset...](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Dataset_generation#about-the-dataset)

-   [The sources](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Dataset_generation#the-sources)
-   [Produce unfiltered VCF with all samples](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Dataset_generation#produce-unfiltered-vcf-with-all-samples)
-   [Filtering of the dataset](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Dataset_generation#filtering-of-the-dataset)

## [Population structure](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Pop_structure#population-structure)

-   [PCA](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Pop_structure#pca)
-   [Ancestry plots](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Pop_structure#ancestry-plots)
-   [Tree of the samples (maximum likelihood)](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Pop_structure#tree-of-the-samples-maximum-likelihood)
-   [IBD Network](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Pop_structure#ibd-network)
-   [Population structure for *P. simium* only](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Pop_structure#population-structure-for-p-simium-only)


## [Population graphs and *f*-statistics](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Pop_graphs_fstats#population-graphs-and-f-statistics)

-   [TreeMix](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Pop_graphs_fstats#treemix)
-   [AdmixtureBayes](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Pop_graphs_fstats#admixturebayes)
-   [*f~4~*-statistics](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Pop_graphs_fstats#f4-statistics)

## [Demographic history](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Demography#demographic-history)

-   [Nucleotide diversity (pi)](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Demography#nucleotide-diversity-pi)
-   [Coalescence with Relate](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Demography#coalescence-with-relate)
-   [Isolation dynamics and migration rate changes over time with MSMC-IM](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Demography#isolation-dynamics-and-migration-rate-changes-over-time-with-msmc-im)

## [Genome scans for selection](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Selection_scan#genome-scans-for-selection)

-   [Haplotype-based tests (*iHS*,*XP-EHH* and *Rsb*)](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Selection_scan#haplotype-based-tests-ihsxp-ehh-and-rsb)
-   [Selection detection with Relate](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Selection_scan#selection-detection-with-relate)
-   [Genome scans with *PBS*](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Selection_scan#genome-scans-with-pbs)
-   [Identification of the genes](https://github.com/MargauxLefebvre/EvoHistory_Simium/tree/main/Selection_scan#identification-of-the-genes)
