import unittest
import os
import subprocess

run_pipeline_tests = True
run_QC_tests = True
run_deseq_tests = True

class Tester(unittest.TestCase):

    # TODO: split this into separate tests for different steps of the pipeline
    def test_pipeline(self):
        if(run_pipeline_tests):
            print('testing pipeline')
            # run snakemake pipeline
            p1 = subprocess.Popen('snakemake -s workflow/Snakefile.temposeq --cores=1 -R all --delete-all-output', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            p1.communicate()
            p2 = subprocess.Popen('snakemake -s workflow/Snakefile.temposeq --cores=1 all', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            p2.communicate()

            # Check exits code and other expected output
            self.assertEqual(0, p1.returncode)
            self.assertEqual(0, p2.returncode)


    def test_QC(self):    # run snakemake pipeline
        if(run_QC_tests):
            print('testing QC')
            # run QC
            p1 = subprocess.Popen('cp data/processed/count_table_mini.csv data/processed/count_table.csv; Rscript scripts/render_studywide_QC_report.R', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            p1.communicate()
            self.assertEqual(0, p1.returncode)

            # check QC output
                  
    def test_deseq(self):
        if(run_deseq_tests):
            print('testing deseq')
            # run DEseq2
            p1 = subprocess.Popen('Rscript scripts/render_DESeq2_report.R', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            stdout, stderr = p1.communicate()
            print(p1.returncode)
            self.assertEqual(0, p1.returncode)

            # check DEseq2 output

if __name__ == '__main__':
    unittest.main()