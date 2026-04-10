Phasing of genotypes for SNPs on the Illumina OmniExpress chip extracted from the GENESTAR whole genome sequence data

by Alexandre Bureau

Markers on the Illumina OmniExpress SNP chip were extracted from bim files for this genotyping array, with physical positions lifted over from hg19 to hg38 using the liftOver function of the R package rtracklayer with hg19ToHg38.over.chain. Variants in the original bim file or in the bim file obtained after flipping strands to conform to the Haplotype Reference Consortium strands using HRC-1000G-check-bim.pl from https://www.well.ox.ac.uk/~wrayner/tools/ were listed in the files 

/dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/OmniExpress_chr*.txt

Variants obtained from the whole genome sequence after the pre-processing steps executed using the script /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/bcftools_plink_pass.sh and described in the /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased/README_pass_SGE_Slurm.txt file are in the files 

/dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/freeze.8.chr*.pass.gtonly.minDP10.dbGaPID.parents.Mendelcompat

1- Variants in these files were extracted if they matched variants in the OmniExpress_chr*.txt files using the script /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/extractOmniExpress.sh using the call

qsub -t 1-22 extractOmniExpress.sh 

Phasing

The phasing per se was performed with Shapeit2 for each autosome. It was performed twice.
A) On the data at this stage of processing
B) On the data after removing variants with likely genotyping errors based on the duoHMM error detection algorithm

Each phasing run involved

i) The inference of the most likely phase indicators. 
ii) The conversion of the phase indicators into nucleotides on each haplotype

2- Initial phasing run

The commands for the phasing steps on the data at this stage of processing are contained in the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased/shapeit_OmniExpress.sh, which was executed for each autosome with the commands

qsub -t 1-22 -pe local 8 -R y -l mem_free=8G,h_vmem=8G shapeit_OmniExpress.sh

3- Identifying genotypes that are unlikely within a duo/trio based on patterns of recombination

This step uses the haps files from the previous step. These haps files need to be uncompressed with bunzip2 before running duoHMM. This is why the maximal file size needs to be raised to 40G.

The commands for this error detection step are contained in the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased/check_geno_error_OmniExpress.sh, which was executed for each autosome with the commands 

qsub -t 1-22 -l h_fsize=40G -l mem_free=31G,h_vmem=32G check_geno_error_OmniExpress.sh

4- Removing variants with likely genotyping errors

Positions where HMM flagged a genotyping error with a probability > 0.95 are extracted, the variant names at that position retrieved from the corresponding bim file, and these variants are then removed from the genotype files using Plink.

The commands for this error cleaning step are contained in the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased/list_variant_errors_OmniExpress.txt, which was executed for each autosome with the commands 

qsub -t 1-22 list_variant_errors_OmniExpress.sh

5- Second phasing run

In this run, graph files containing the graphical representation of the haplotype structure of the sample is written in uncompressed files. 

The commands for the phasing steps on the cleaned data are contained in the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased/shapeit_clean.sh, which was executed for each autosome with the command

qsub -t 1-22 -pe local 8 -R y -l mem_free=8G,h_vmem=8G shapeit_OmniExpress_clean.sh

6- Recombination inference

Ten haplotype sets are generated from the graph with different seeds using Shapeit2, the recombination map of each meiosis in the haplotype sets is calculated using duoHMM and the maps are averaged.

The commands for this recombination inference are contained in the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/recombinations/recombination_inference_OmniExpress_clean.sh, which was executed for each autosome with the command

qsub -t 1-22 -l mem_free=16G,h_vmem=16G recombination_inference_OmniExpress_clean.sh

The final output is contained in the files 

/dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/recombinations/chr*.OmniExpress.clean.recombinations
