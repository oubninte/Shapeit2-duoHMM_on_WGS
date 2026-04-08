Phasing of the GENESTAR whole genome sequence SNVs with pass status

Running scripts were updated to the Slurm scheduler to add chromosome 9 in 2026.

by Alexandre Bureau

The phasing was performed on the bcf files freeze.8.chr*.pass_and_fail.gtonly.minDP10.bcf


The initial processing involved for each autosome

1- Removing monomorphic sites from the bcf files and keeping only variants with FILTER=PASS using bcftools

2- Converting bcf files into bed/bim/fam files, retaining only single nucleotide variants (SNVs) while applying QC steps based on Hardy-Weinberg equilibrium test and proportion of missing values (5%) using Plink

3- Recoding fam files with Plink, including updating the family IDs of the subjects
This step requires the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/WGSupdateFID.txt

4- Adding the father id, mother id and sex of the subjects to the fam file with Plink
This step requires the files  /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/WGSupdateParents.txt and /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/WGSupdateSex.txt

5- Checking for and correcting Mendelian inconsistencies using Plink, and assigning ids to variants in the format chr:position:A1:A2


The commands for these processing steps are contained in the script /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/bcftools_plink_pass.sh, which can be executed for the 22 autosomes with the call

qsub -t 1-22 bcftools_plink_pass.sh

6- Trimming of markers outside the genetic map

To solve errors of the type

ERROR BAD RECOMBINATION RATE AT MARKER 4049464 rho = nan
223.457 -nan

markers outside the genetic map were removed.

The command for this step is contained in the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/maprange_pass.sh. This was performed for chromosomes 1 to 4, 7, 11, 17 and 18. Chromosome 9 was split into its p and q arms to avoid failing the initial phasing run. The p arm was renamed chromosome 29. Phasing with Shapeit2 was successful for the other autosomes without trimming. The data files were overwritten when performing this step.

qsub -t 1-4 maprange_pass.sh 
qsub -t 7 maprange_pass.sh 
sbatch --array=9 maprange_pass_slurm.sh 
sbatch --array=29 maprange_pass_slurm.sh 
qsub -t 11 maprange_pass.sh 
qsub -t 17-18 maprange_pass.sh 


Phasing

The phasing per se was performed with Shapeit2 for each autosome. It was performed twice.
A) On the data at this stage of processing
B) On the data after removing variants with likely genotyping errors based on the duoHMM error detection algorithm

Each phasing run involved

i) The inference of the most likely phase indicators. This process is very memory consuming.
ii) The conversion of the phase indicators into nucleotides on each haplotype

8- Initial phasing run

The commands for the phasing steps on the data at this stage of processing are contained in the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased/shapeit_pass.sh, which was executed for each autosome with the commands

qsub -t 1-2 -l h_fsize=40G -pe local 8 -R y -l mem_free=16G,h_vmem=16G shapeit_pass.sh
qsub -t 3-22 -l h_fsize=40G -pe local 8 -R y -l mem_free=8G,h_vmem=8G shapeit_pass.sh

Chromosome 9 failed this step. The solution was to split chromosome 9 into its p and q arms. The p arm was renamed chromosome 29. This required to copy the map file dcl01/mathias1/data/TOPMED_WGS/GRCh38_genetic_maps/shapeit.chr9.GRCh38.map as dcl01/mathias1/data/TOPMED_WGS/GRCh38_genetic_maps/shapeit.chr29.GRCh38.map

sbatch --array=9 --mem-per-cpu=5G --cpus-per-task=8 --time=4-00:00 shapeit_pass_slurm.sh
sbatch --array=29 --mem-per-cpu=5G --cpus-per-task=8 --time=2-00:00 shapeit_pass_slurm.sh


9- Identifying genotypes that are unlikely within a duo/trio based on patterns of recombination

This step uses the haps files from the previous step. These haps files need to be uncompressed with bunzip2 before running duoHMM. This is why the maximal file size needs to be raised to 40G.

The commands for this error detection step are contained in the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased/check_geno_error_pass.sh, which was executed for each autosome with the commands 

qsub -t 1-5 -l h_fsize=40G -l mem_free=24G,h_vmem=24G check_geno_error_pass.sh
qsub -t 6-8 -l h_fsize=40G -l mem_free=16G,h_vmem=16G check_geno_error_pass.sh
sbatch --array=9 --mem-per-cpu=4G --cpus-per-task=4 check_geno_error_pass_slurm.sh 
sbatch --array=29 --mem-per-cpu=4G --cpus-per-task=4 check_geno_error_pass_slurm.sh 
qsub -t 10-17 -l h_fsize=40G -l mem_free=16G,h_vmem=16G check_geno_error_pass.sh
qsub -t 18-22 -l h_fsize=40G -l mem_free=8G,h_vmem=8G check_geno_error_pass.sh


10- Removing variants with likely genotyping errors

Positions where HMM flagged a genotyping error with a probability > 0.95 are extracted, the variant names at that position retrieved from the corresponding bim file, and these variants are then removed from the genotype files using Plink.

The commands for this error cleaning step are contained in the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased/list_variant_errors_pass.txt, which was executed for each autosome with the commands 

qsub -t 1-8 list_variant_errors_pass.txt
sbatch --array=9 list_variant_errors_pass_Slurm.txt
sbatch --array=29 list_variant_errors_pass_Slurm.txt
qsub -t 10-22 list_variant_errors_pass.txt

11- Second phasing run

In this run, graph files containing the graphical representation of the haplotype structure of the sample is written in uncompressed files. This is why the maximal file size needs to be raised to 40G.

The commands for the phasing steps on the cleaned data are contained in the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased/shapeit_pass_clean.sh, which was executed for each autosome with the commands

qsub -t 1-2 -l h_fsize=40G -pe local 8 -R y -l mem_free=16G,h_vmem=16G shapeit_pass_clean.sh
qsub -t 3-8 -l h_fsize=40G -pe local 8 -R y -l mem_free=8G,h_vmem=8G shapeit_pass_clean.sh
sbatch --array=9 --mem-per-cpu=5G --cpus-per-task=8 --time=4-00:00 shapeit_pass_clean_Slurm.sh
sbatch --array=29 --mem-per-cpu=5G --cpus-per-task=8 --time=4-00:00 shapeit_pass_clean_Slurm.sh
qsub -t 10-22 -l h_fsize=40G -pe local 8 -R y -l mem_free=8G,h_vmem=8G shapeit_pass_clean.sh

12- Recombination inference

Ten haplotype sets are generated from the graph with different seeds using Shapeit2, the recombination map of each meiosis in the haplotype sets is calculated using duoHMM and the maps are averaged.

The commands for this recombination inference are contained in the file /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/recombinations/recombination_inference_pass_clean.sh, which was executed for each autosome with the commands

qsub -t 1-6 -l mem_free=31G,h_vmem=32G recombination_inference_pass_clean.sh
qsub -t 7-8 -l mem_free=16G,h_vmem=16G recombination_inference_pass_clean.sh
sbatch --array=9 --mem-per-cpu=4G --cpus-per-task=4 --time=4-00:00 recombination_inference_pass_clean_Slurm.sh 
sbatch --array=29 --mem-per-cpu=4G --cpus-per-task=4 --time=4-00:00 recombination_inference_pass_clean_Slurm.sh 
qsub -t 10-17 -l mem_free=16G,h_vmem=16G recombination_inference_pass_clean.sh
qsub -t 18-22 -l mem_free=8G,h_vmem=8G recombination_inference_pass_clean.sh

The phased haplotype files are /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased/freeze.8.chr*.pass.gtonly.minDP10.dbGaPID.parents.Mendelcompat.clean.phased.hap.bz2

The final output is contained in the files 

/dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/recombinations/chr*.pass.clean.recombinations
