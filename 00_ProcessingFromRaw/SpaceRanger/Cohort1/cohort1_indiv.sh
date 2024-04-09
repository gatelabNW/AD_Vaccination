#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition normal
#SBATCH --job-name spaceranger
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 128G
#SBATCH --time 20:00:00
#SBATCH --output /projects/p31535/alex/AN1792//AN1792_cohort1_distance/logs/%x.oe%j.log
#SBATCH --verbose

# Specify arguments
# $1 = sample_id
# $2 = sample_name
# $3 = slide serial number
# $4 = slide area

# Define directories
alex_dir="/projects/b1042/Gate_Lab/alex/submit-spaceranger/AN1792_cohort1_distance"
thomas_dir="/projects/p31535/thomas/SpaceRanger/spaceranger-2.1.1/"
ref_dir="/projects/p31535/spatial-transcriptomics/spacerange-ref-files-before-transfer"
cytassist_image_dir="${alex_dir}/cytassist/"
hi_res_image_dir="${alex_dir}/high-res/"
alignment_dir="${alex_dir}/json/"
fastq_dir="${alex_dir}/fastq/${2}/"
output_dir="/projects/b1042/Gate_Lab/alex/spaceranger-output/AN1792/AN1792_cohort1_distance/"

# Navigate to output directory
cd $output_dir

# Run spaceranger
/projects/p31535/thomas/SpaceRanger/spaceranger-2.1.1/spaceranger count \
--id="$1" \
--description="AN1792_cohort1_distance" \
--transcriptome="/projects/p31535/spatial-transcriptomics/spacerange-ref-files/refdata-gex-GRCh38-2020-A" \
--probe-set="${thomas_dir}probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv" \
--fastqs="${fastq_dir}" \
--cytaimage="${cytassist_image_dir}$1.tif" \
--darkimage="${hi_res_image_dir}${2}_distance_multipage.tif" \
--loupe-alignment="${alignment_dir}${3}-${4}.json" \
--slide="${3}" \
--area="${4}" \
--localcores=16 \
--localmem=128
