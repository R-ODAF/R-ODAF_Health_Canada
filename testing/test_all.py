import os
import subprocess
import pytest

run_pipeline_tests = True
run_QC_tests = True
run_deseq_tests = True


@pytest.mark.dependency()
def test_pipeline():
    if(run_pipeline_tests):
        print('testing pipeline')
        # run snakemake pipeline
        p1 = subprocess.Popen('snakemake -s workflow/Snakefile.temposeq --cores=1 -R all --delete-all-output', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        p1.communicate()
        p2 = subprocess.Popen('snakemake -s workflow/Snakefile.temposeq --cores=1 all', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        p2.communicate()

        # Check exits code and other expected output
        assert p1.returncode == 0
        assert p2.returncode == 0


@pytest.mark.dependency(depends=["test_pipeline"])
def test_QC():    # run snakemake pipeline
    if(run_QC_tests):
        print('testing QC')
        # run QC
        p1 = subprocess.Popen('cp data/processed/count_table_mini.csv data/processed/count_table.csv; Rscript scripts/render_studywide_QC_report.R', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        p1.communicate()
        assert p1.returncode == 0

        # check QC output

@pytest.mark.dependency(depends=["test_QC"])
def test_deseq():
    if(run_deseq_tests):
        print('testing deseq')
        # run DEseq2
        p1 = subprocess.Popen('Rscript scripts/render_DESeq2_report.R', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        p1.communicate()
        assert p1.returncode == 0

        # check DEseq2 output
