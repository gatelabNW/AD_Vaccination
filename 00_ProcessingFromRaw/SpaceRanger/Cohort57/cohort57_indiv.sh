#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name spaceranger
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 128G
#SBATCH --time 20:00:00
#SBATCH --output /projects/p31535/alex/AN1792_cohort57_plaques_vessels/logs/%x.oe%j.log
#SBATCH --verbose

# Specify arguments
# $1 = csv
# $2 = sample_name
# $3 = slide serial number
# $4 = slide area

# Define directories
alex_dir="/projects/b1042/Gate_Lab/alex/submit-spaceranger/AN1792_vesselDist_plaqueDist"
thomas_dir="/projects/p31535/thomas/SpaceRanger/spaceranger-2.1.1/"
ref_dir="/projects/p31535/spatial-transcriptomics/spacerange-ref-files"
cytassist_image_dir="${alex_dir}/cytassist/"
hi_res_image_dir="${alex_dir}/high-res/"
alignment_dir="${alex_dir}/json/"
fastq_dir="${alex_dir}/fastq/${2}/"
output_dir="/projects/b1042/Gate_Lab/alex/spaceranger-output/AN1792/AN1792_vesselDist_plaqueDist/"

# Navigate to output directory
cd $output_dir

# Run spaceranger
/projects/p31535/thomas/SpaceRanger/spaceranger-2.1.1/spaceranger count \
--id="$2" \
--description="AN1792_cohort57_with_vesselDM_and_plaqueDM" \
--transcriptome="/projects/p31535/spatial-transcriptomics/spacerange-ref-files/refdata-gex-GRCh38-2020-A" \
--probe-set="${thomas_dir}probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv" \
--feature-ref="${ref_dir}/spaceranger-2.1.1/feature_refs/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv" \
--libraries="${alex_dir}/csv/${1}_IF.csv" \
--cytaimage="${cytassist_image_dir}$1.tif" \
--darkimage="${hi_res_image_dir}${2}_multipage_orig-dapi-vesselDM-plaqueDM-rawBleachedIba1.tif" \
--loupe-alignment="${alignment_dir}${3}-${4}.json" \
--slide="${3}" \
--area="${4}" \
--localcores=16 \
--localmem=128
