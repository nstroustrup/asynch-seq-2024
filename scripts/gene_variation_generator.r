
das = data.frame(Name = c("N2_ZC443.3_Set5_Day8",
"N2_D1086.7_Set5_Day8",
"N2_atp5_Set5_Day8",
"N2_nlp31_Set5_Day8",
"N2_col122_Set3_Day8",
"N2_nlp28_Set3_Day8",
"N2_atp5_Set5_Day8",
"N2_nlp15_Set3_Day8",
"N2_col122_Set3_Day8",
"N2_C25E10.8_Set2_Day8",
"N2_pat10_Set3_Day8"),
gene_name=c("ZC443.3","D1086.7","atp-5","nlp-31","col-122","nlp-28","atp-5","nlp-15","col-122","C25E10.8","pat-10"),
hit_subject = c(rep("d1_sw",5),rep("d8_sw",6)))


output_header = r"(#!/bin/bash
#$ -l virtual_free=32G,h_rt=6:00:00
#$ -N ns_sw_FF   #job name
#$ -q short-sl7  #queue
#$ -t 1-50
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/gene_variability_single_worms.R -n 250 -i $SGE_TASK_ID -r F -c F --split_bootstraps_across_nodes 50 --differential_analysis_comparing_da_model T)";

output_header_step_2 = r"(#!/bin/bash
#$ -l virtual_free=64G,h_rt=24:00:00
#$ -N ns_sw_TF   #job name
#$ -q long-sl7  #queue
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/gene_variability_single_worms.R -n 250 -r T -c F --split_bootstraps_across_nodes 50 --differential_analysis_comparing_da_model T)";

k = 1;
for (model_name in unique(das$hit_subject)){
	for (i in 1:nrow(das)){
	da_name = das$Name[i]
	da_label = str_replace(das$gene_name[i],"-","")
	gene_name = das$gene_name[i]
	job_suffix = paste("-m",model_name,
			paste0("--differential_analysis_model \"data/differential_analysis/tables/net_validation/",da_name,".csv.gz\""),
			"--differential_analysis_model_name",paste0("\"",da_label,"\""),
			paste0("--differential_analysis_model_foldchange_column \"RNAi_",da_label,"_vs_EV_log2FoldChange\""),
			"--differential_analysis_RNAi_GeneSymbol_to_remove",paste0("\"",gene_name,"\"")
			)
	
	script = paste(output_header,job_suffix);
	cat(script,file=paste0("scripts/gene_variability_scripts/run_",k,".sh"))
	
	script_step_2 = paste(output_header_step_2,job_suffix);
	cat(script_step_2,file=paste0("scripts/gene_variability_scripts/run_step2_",k,".sh"))
	k=k+1;
	#stop()
	}
}
		