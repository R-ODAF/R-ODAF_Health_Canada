import os

run_pipeline = True
run_QC = True
run_deseq = True

cwd = os.getcwd()
print(cwd)

os.chdir("..") # should be main dir

cwd = os.getcwd()
print(cwd)


# run snakemake pipeline
os.system('snakemake --cores=1 -s workflow/Snakefile.temposeq -R all_temposeq')

# check pipeline output


# run QC
os.system('Rscript scripts/render_studywide_QC_report.R')

# check QC output


# run DEseq2
os.system('Rscript scripts/render_DESeq2_report.R')

# check DEseq2 output