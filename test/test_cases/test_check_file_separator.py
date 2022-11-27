import unittest
import subprocess as sp
from test.utils import run_snakemake

class Test_check_file_separator(unittest.TestCase):
    ''' '''

    def test_tsv_file(self):
        p = run_snakemake('snakemake -p -c 1 results/check-file-separator_checked_separator.tsv.gz --force')

        stdout, stderr = p.communicate()

        with self.subTest(msg = 'Return code check'):
            self.assertEqual(0, p.returncode)

    def test_csv_file(self):
        p = run_snakemake('snakemake -p -c 1 results/check-file-separator_checked_separator.csv.gz --force')

        stdout, stderr = p.communicate()

        with self.subTest(msg = 'Return code check'):
            self.assertEqual(0, p.returncode)

    def test_semicolon_file(self):
        p = run_snakemake('snakemake -p -c 1 results/check-file-separator_checked_separator.semicolon.gz --force')

        stdout, stderr = p.communicate()

        with self.subTest(msg = 'Return code check'):
            self.assertEqual(0, p.returncode)

    def test_space_file(self):
        p = run_snakemake('snakemake -p -c 1 results/check-file-separator_checked_separator.space.gz --force')

        stdout, stderr = p.communicate()

        with self.subTest(msg = 'Return code check'):
            self.assertEqual(0, p.returncode)

if __name__ == '__main__':
    unittest.main()
