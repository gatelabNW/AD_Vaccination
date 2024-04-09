#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name spaceranger
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 128G
#SBATCH --time 4:00:00
#SBATCH --output /projects/p31535/alex/AN1792/logs/%x_oe%j.log
#SBATCH --verbose

# Define directories
main_dir="/projects/p31535/spatial-transcriptomics"
alex_dir="/projects/b1169/alex"
ref_dir="${main_dir}/spacerange-ref-files/"
cytassist_image_dir="${alex_dir}/submit-spaceranger/AN1792_cytassist_images/"
hi_res_image_dir="${alex_dir}/submit-spaceranger/Full_DAB_H_amyloid/AN1792_binary_with_amyloid_grayscale_multipage/8-bit/"
alignment_dir="${alex_dir}/submit-spaceranger/Full_DAB_H_amyloid/AN1792_binary_with_amyloid_grayscale_jsons/8-bit/"
fastq_dir="/projects/b1169/alex/submit-spaceranger/AN1792_fastqs/102-1/"
output_dir="/projects/b1169/alex/spaceranger-output/AN1792/Binary/"

# Navigate to output directory
cd $output_dir

# Run spaceranger
/projects/p31535/spatial-transcriptomics/spacerange-ref-files/spaceranger-2.1.1/spaceranger count \
--id="AN1792_102-1" \
--sample="AN1792_102-1" \
--description="AN1792_102-1_binary" \
--transcriptome="${ref_dir}refdata-gex-GRCh38-2020-A" \
--probe-set="${ref_dir}Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv" \
--fastqs="${fastq_dir}" \
--cytaimage="${cytassist_image_dir}AN1792 102-1.tif" \
--colorizedimage="${hi_res_image_dir}102-1_multi_SR_8bit.tif" \
--loupe-alignment="${alignment_dir}V42Y23-312-A1.json" \
--slide="V42Y23-312" \
--area="A1" \
--localcores=16 \
--localmem=128