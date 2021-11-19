import unittest
import os
import subprocess

run_pipeline = True
run_QC = False
run_deseq = False

# cwd = os.getcwd()
# print(cwd)

# os.chdir("..") # should be main dir

# cwd = os.getcwd()
# print(cwd)

class Tester(unittest.TestCase):

    def test_pipeline(self):
        if(run_pipeline):
            # run snakemake pipeline
            p = subprocess.Popen('snakemake -s workflow/Snakefile.temposeq --cores=1 -R all --delete-all-output', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            p.communicate()
            p = subprocess.Popen('snakemake -s workflow/Snakefile.temposeq --cores=1 all', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            p.communicate()
            
            # Check exits code and other expected output            
            self.assertEqual(0, p.returncode)
            #self.assertTrue(os.path.isfile('some.pipeline.output.txt'))


    def test_QC(self):    # run snakemake pipeline
        if(run_deseq):
            # run QC
            os.system('Rscript scripts/render_studywide_QC_report.R')

            # check QC output

    def test_deseq(self):
        if(run_deseq):
            # run DEseq2
            os.system('Rscript scripts/render_DESeq2_report.R')

            # check DEseq2 output

if __name__ == '__main__':
    unittest.main()