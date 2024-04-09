#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition normal
#SBATCH --job-name spaceranger
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 128G
#SBATCH --time 20:00:00
#SBATCH --output /projects/p31535/alex/AN1792_IF/logs/%x.oe%j.log
#SBATCH --verbose

# Specify arguments
# $1 = sample_name
# $2 = sample_id
# $3 = slide serial number
# $4 = slide area
# $5 = high res iba1
# $6 = high res amyloid
# $7 = cytassist
# $8 = csv

# Define directories
alex_dir="/projects/b1042/Gate_Lab/alex/submit-spaceranger/AN1792_cohort3"
thomas_dir="/projects/p31535/thomas/SpaceRanger/spaceranger-2.1.1/"
ref_dir="/projects/p31535/spatial-transcriptomics/spacerange-ref-files"
cytassist_image_dir="${alex_dir}/cytassist/"
hi_res_image_dir="${alex_dir}/high-res/iba1_w_rollingball/"
alignment_dir="${alex_dir}/json/iba1-redone/"
# fastq_dir="${alex_dir}/fastq/${2}/"
output_dir="/projects/b1042/Gate_Lab/alex/spaceranger-output/AN1792/cohort3/iba1_w_rollingball/"

# Navigate to output directory
cd $output_dir

# Run spaceranger
/projects/p31535/thomas/SpaceRanger/spaceranger-2.1.1/spaceranger count \
--id="$1" \
--description="AN1792_cohort3" \
--transcriptome="/projects/p31535/spatial-transcriptomics/spacerange-ref-files/refdata-gex-GRCh38-2020-A" \
--probe-set="${thomas_dir}probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv" \
--feature-ref="${ref_dir}/spaceranger-2.1.1/feature_refs/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv" \
--libraries="${alex_dir}/csv/$8.csv" \
--cytaimage="${cytassist_image_dir}$7.tiff" \
--darkimage="${hi_res_image_dir}$5_RB.tif" \
--loupe-alignment="${alignment_dir}${3}-${4}.json" \
--slide="${3}" \
--area="${4}" \
--localcores=16 \
--localmem=128
