#!/bin/bash
#SBATCH -J template
#SBATCH -o out
#SBATCH -e err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=80G         # memory per cpu-core (4G is default)
#SBATCH --time 1-00:00:00        # DAYS-HOURS:MINUTES:SECONDS

source /home/bjarnold/miniforge3/etc/profile.d/conda.sh
conda activate cyvcf

# VCF="/scratch/gpfs/bjarnold/chernobylWolves/data/VCFs/01_depth_missdata_filtered/FOR_LUDWIG_DEADLINE/01_merged_vcf_new_samples_with_outgroup/new_samples_with_outgroup_filtered.vcf.gz"
#VCF="/scratch/gpfs/bjarnold/chernobylWolves/data/VCFs/01_depth_missdata_filtered/FOR_LUDWIG_DEADLINE/01_merged_vcf_new_samples_with_outgroup/test.vcf"
#VCF="/scratch/gpfs/bjarnold/chernobylWolves/data/VCFs/01_depth_missdata_filtered/DP8/concat.vcf.gz"

OUT_DIR="full_vcf"
mkdir -p ${OUT_DIR}

POP_FILE_DIR="/scratch/gpfs/bjarnold/chernobylWolves/analyses/sliding_window/pop_files"
VCF="/scratch/gpfs/cl5291/GW_newDir/14.SNPeff/snpEff/WG_wolf_concat_DP8_annotated.vcf"
# VCF="./test.vcf"
WIN=21
OFF=1
MINAF=1

python daf_fst_faywuh.py \
--vcf ${VCF} \
--winsize ${WIN} \
--offset ${OFF} \
--min_AF ${MINAF} \
--pop_file ${POP_FILE_DIR}"/CEZ_BLR_YLS.txt" \
--focal_pop "CEZ" \
--sister_pop "BLR" \
--outgroup_pop "YLS" \
--output_fst ${OUT_DIR}/"CEZ_BLR_YLS_WIN${WIN}_OFF${OFF}_MINAF${MINAF}_FST.txt" \
--output_faywuh ${OUT_DIR}/"CEZ_BLR_YLS_WIN${WIN}_OFF${OFF}_MINAF${MINAF}_FayWuH.txt" \
--output_ddaf ${OUT_DIR}/"CEZ_BLR_YLS_WIN${WIN}_OFF${OFF}_MINAF${MINAF}_DDAF.txt"


