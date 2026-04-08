module load bcftools

bcftools view /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/minDP10/freeze.8.chr${SGE_TASK_ID}.pass_and_fail.gtonly.minDP10.bcf -f .,PASS -c 1:minor -Ob > /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/chr${SGE_TASK_ID}.pass.bcf

/users/abureau/plink_linux_x86_64_20200103/plink --bcf /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/chr${SGE_TASK_ID}.pass.bcf --allow-extra-chr --snps-only 'just-acgt' --remove /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/WGStoremove.txt --hwe 1e-06 --geno 0.05 --make-bed --out /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/freeze.8.chr${SGE_TASK_ID}.pass.gtonly.minDP10

/users/abureau/plink_linux_x86_64_20200103/plink --bfile /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/freeze.8.chr${SGE_TASK_ID}.pass.gtonly.minDP10 --list-duplicate-vars --update-ids /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/WGSupdateFID.txt --make-bed --out /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/freeze.8.chr${SGE_TASK_ID}.pass.gtonly.minDP10.dbGaPID

/users/abureau/plink_linux_x86_64_20200103/plink --bfile /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/freeze.8.chr${SGE_TASK_ID}.pass.gtonly.minDP10.dbGaPID --update-parents /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/WGSupdateParents.txt --update-sex /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/WGSupdateSex.txt --make-bed --out /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/freeze.8.chr${SGE_TASK_ID}.pass.gtonly.minDP10.dbGaPID.parents

/users/abureau/plink_linux_x86_64_20200103/plink --bfile /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/freeze.8.chr${SGE_TASK_ID}.pass.gtonly.minDP10.dbGaPID.parents --set-missing-var-ids @:#:\$1:\$2 --set-me-missing --mendel-duos --mendel-multigen --make-bed --out /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/freeze.8.chr${SGE_TASK_ID}.pass.gtonly.minDP10.dbGaPID.parents.Mendelcompat
