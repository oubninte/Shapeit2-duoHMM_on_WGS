# Extract position of variants with probable genotyping error
perl -ne '@vec = split; print "$vec[3]:$vec[4]\n" if $vec[5]>0.95' < /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased_OmniExpress/freeze.8.chr${SGE_TASK_ID}.pass_and_fail.gtonly.minDP10.dbGaPID.parents.Mendelcompat.OmniExpress.phased.generr | sort | uniq > /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/chr${SGE_TASK_ID}.exclude 

# Extract IDs (position and alleles) of all variants on the chromosome
awk '{print $2}' /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/freeze.8.chr${SGE_TASK_ID}.pass_and_fail.gtonly.minDP10.dbGaPID.parents.Mendelcompat.OmniExpress.bim > /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/chr${SGE_TASK_ID}.ids

# Extract the variants at the positions to exclude
snps=($(awk '{print $1}' /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/chr${SGE_TASK_ID}.exclude))
for isnp in "^${snps[@]}"
do
grep ${isnp} /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/chr${SGE_TASK_ID}.ids >> /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/chr${SGE_TASK_ID}.exclude.ids
done

# Exclude variants with probable genotyping error from Plink files
/users/abureau/plink_linux_x86_64_20200103/plink --bfile /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/freeze.8.chr${SGE_TASK_ID}.pass_and_fail.gtonly.minDP10.dbGaPID.parents.Mendelcompat.OmniExpress --exclude /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/chr${SGE_TASK_ID}.exclude.ids --make-bed --out /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/freeze.8.chr${SGE_TASK_ID}.pass_and_fail.gtonly.minDP10.dbGaPID.parents.Mendelcompat.OmniExpress.clean

