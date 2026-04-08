for i in {1..10};
do
  /users/abureau/shapeit/bin/shapeit -convert --input-graph /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/family_phased/freeze.8.chr${SGE_TASK_ID}.pass.gtonly.minDP10.dbGaPID.parents.Mendelcompat.clean.graph --output-sample /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/recombinations/chr${SGE_TASK_ID}.pass.clean.sim${i} --include-ind /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/recombinations/fam10seq.inds --seed ${i}
  /users/abureau/shapeit/bin/duohmm -H /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/recombinations/chr${SGE_TASK_ID}.pass.clean.sim${i} -M /dcl01/mathias1/data/TOPMED_WGS/GRCh38_genetic_maps/shapeit.chr${SGE_TASK_ID}.GRCh38.map -R /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/recombinations/chr${SGE_TASK_ID}.pass.clean.sim${i}.rec
done

python /users/abureau/shapeit/bin/mapavg.py /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/recombinations/chr${SGE_TASK_ID}.pass.clean.*.rec > /dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/recombinations/chr${SGE_TASK_ID}.pass.clean.recombinations