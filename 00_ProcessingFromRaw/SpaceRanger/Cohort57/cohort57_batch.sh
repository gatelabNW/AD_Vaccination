# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                      AN1792-Vacc Project                           -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-18-2023
# Written by: Natalie Piehl
# Used on: 12-5-23
# Used by: Alex Edwards
# Summary: Submit spaceranger job on all samples
#
# ------------------------------------------------------------------------------

# Print date
date

# Define directories
main_dir="/projects/b1042/Gate_Lab/alex/"
metadata="${main_dir}submit-spaceranger/AN1792_cohort57_bleached_iba1/AN1792_cohort57_vessels_plaques.csv"
progress="Found metadata"

echo $progress

# Create array of sample names from metadata
sample_arr=()
while IFS=',' read -ra row; do
  sample_arr+=("${row[0]}")
done < $metadata
progress="Created array"

echo $progress

# Remove header names
sample_arr=("${sample_arr[@]:2}")

# For each sample...
counter=0
echo $"counter: $counter"
# for sample in "${sample_arr[0]}"
for sample in "${sample_arr[@]}"
do
    # Update counter and print status
    let counter++
    echo "Processing $counter: $sample now"

    # Get output sample name, slide ID, and slide slot from metadata
    csv=$(awk -v sam="$sample" -F, '{ if ($1 == sam) print $1}' "$metadata")
    echo $csv
    sample_name=$(awk -v sam="$sample" -F, '{ if ($1 == sam) print $2}' "$metadata")
    echo $sample_name
    slide_id=$(awk -v sam="$sample" -F, '{ if ($1 == sam) print $3}' "$metadata")
    echo $slide_id
    slide_slot=$(awk -v sam="$sample" -F, '{ if ($1 == sam) print $4}' "$metadata")
    echo $slide_slot

    # Run spaceranger
    sbatch /projects/b1169/alex/scripts/AN1792/AN1792_IF/cohort57_indiv.sh \
     "$csv" "$sample_name" "$slide_id" "$slide_slot"

    echo "Submitted $counter: $sample"

    # Sleep to give scheduler a break
    sleep 3
done