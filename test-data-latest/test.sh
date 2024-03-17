#!/bin/bash
# Provide one argument - the branch to test. If no branch is provided, it will use main.
if [ $# -eq 0 ]; then
    branch="main"
fi

if [ $# -gt 1 ]; then
    echo "Please provide 1 argument, the branch to checkout."
    exit 1
fi

if [ $# -eq 1 ]; then
    branch="$1"
fi

git clone https://github.com/R-ODAF/R-ODAF_Health_Canada.git

cd R-ODAF_Health_Canada
git checkout ${branch}
git clone https://github.com/EHSRB-BSRSE-Bioinformatics/test-data
rm -r data # Remove existing directory before replacing w/ test data
rm -r config # Remove existing directory before replacing w/ test data
mv test-data/temposeq/* ./
wget https://github.com/EHSRB-BSRSE-Bioinformatics/unify_temposeq_manifests/raw/main/output_manifests/Human_S1500_1.2_standardized.csv

# Run this if you don't already have an environment to use...
#conda env create -f environment.yml
#mamba activate R-ODAF
eval "$(conda shell.bash hook)"
conda activate snakemake
snakemake --cores 8 --use-conda
Rscript scripts/render_studywide_QC_report.R
Rscript scripts/run_DESeq2.R
Rscript scripts/render_DESeq2_report.R

md5sum data/processed/count_table.tsv > checksum.count_table
md5sum data/metadata/metadata.QC_applied.txt > checksum.metadata.QC_applied
md5sum analysis/QC/MultiQC_Report_data.zip > checksum.multiQC
md5sum analysis/QC/details/QC_per_sample.txt > checksum.QC_per_sample
md5sum analysis/DEG_lists/BaP/*.txt > checksums.BaP_DEGs
md5sum analysis/DEG_lists/CISP/*.txt > checksums.CISP_DEGs
md5sum analysis/pathway_analysis/BaP/*WikiPathways* > checksums.BaP_pathways
md5sum analysis/pathway_analysis/CISP/* > checksums.CISP_pathways
md5sum analysis/BMD_and_biomarker_files/*.txt > checksums.BMD_files
md5sum analysis/BMD_and_biomarker_files/*/*.txt >> checksums.BMD_files

cmp truth_checksums/checksum.count_table checksum.count_table
cmp truth_checksums/checksum.metadata.QC_applied checksum.metadata.QC_applied
cmp truth_checksums/checksum.QC_per_sample checksum.QC_per_sample
echo "BaP DEG checksums from this run:"
cat checksums.BaP_DEGs # Can be informative if there are errors
echo "BaP DEG checksums from 'truth' set:"
cat truth_checksums/checksums.BaP_DEGs # Can be informative if there are errors
echo "CISP DEG checksums from this run:"
cat checksums.CISP_DEGs
echo "CISP DEG checksums from 'truth' set:"
cat truth_checksums/checksums.CISP_DEGs # Can be informative if there are errors
echo "BaP pathway checksums from this run:"
cat checksums.BaP_pathways
echo "BaP pathway checksums from 'truth' set:"
cat truth_checksums/checksums.BaP_pathways

diff -q <(sort -u checksums.BaP_DEGs | awk {'print $1'}) \
<(grep -Fxf \
<(cat checksums.BaP_DEGs | awk {'print $1'}) \
<(cat truth_checksums/checksums.BaP_DEGs | awk {'print $1'} | sort -u))

diff -q <(sort -u checksums.CISP_DEGs | awk {'print $1'}) \
<(grep -Fxf \
<(cat checksums.CISP_DEGs | awk {'print $1'}) \
<(cat truth_checksums/checksums.CISP_DEGs | awk {'print $1'} | sort -u))

cmp truth_checksums/checksums.BaP_pathways checksums.BaP_pathways
cat checksums.BMD_files # Can be informative if there are errors
cat truth_checksums/checksums.BMD_files # Can be informative if there are errors
cmp truth_checksums/checksums.BMD_files checksums.BMD_files

