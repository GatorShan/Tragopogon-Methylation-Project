import multiprocessing
import subprocess
import shlex
from methylpy.call_mc_pe import merge_sorted_multimap_pe, call_methylated_sites_pe # add _pe
from methylpy.call_mc_se import remove_clonal_bam
import os

num_procs = 16
reference_fasta = "/ufrc/soltis/shan158538/Methylation/OutPut/Build_reference/Tdub_lambda.fasta"
sample = "T.dubius_3040-6-2" # CHECK
path_to_output = "/orange/soltis/shan158538/methylpy_output/methylpy_dataprocessing_Tdu_3040-6-2_V1/" # CHECK

# Calling methylated sites
output_bam_file = path_to_output + sample + "_libA_processed_reads_no_clonal.bam"
call_methylated_sites_pe(output_bam_file,
                         sample,
                         reference_fasta,
                         unmethylated_control="NC_001416.1_lambda:",
                         sig_cutoff=0.01,
                         num_procs=num_procs,
                         num_upstr_bases=0,
                         num_downstr_bases=2,
                         generate_mpileup_file=True,
                         compress_output=True,
                         bgzip=False,
                         path_to_bgzip="",
                         path_to_tabix="",
                         min_cov=2,
                         binom_test=True,
                         remove_chr_prefix=True,
                         sort_mem="500M",
                         path_to_files=path_to_output,
                         path_to_samtools="",
                         min_base_quality=1,
                         keep_temp_files=True)


