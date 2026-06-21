import unittest
import os
import sys
import tempfile
import gzip

# Add root directory to path to allow import
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

class TestPISEModules(unittest.TestCase):
    def test_filter_reads_imports(self):
        try:
            from pise import filter_reads
            self.assertTrue(hasattr(filter_reads, 'load_primers'))
            self.assertTrue(hasattr(filter_reads, 'process_forward_reads'))
            self.assertTrue(hasattr(filter_reads, 'process_batch_worker'))
        except ImportError as e:
            self.fail(f"Could not import pise.filter_reads: {e}")

    def test_qc_imports(self):
        try:
            from pise import qc
            self.assertTrue(hasattr(qc, 'run_qc'))
            self.assertTrue(hasattr(qc, 'load_known_indices'))
            self.assertTrue(hasattr(qc, 'get_or_create_handle'))
        except ImportError as e:
            self.fail(f"Could not import pise.qc: {e}")

    def test_extract_pairs_imports(self):
        try:
            from pise import extract_pairs
            self.assertTrue(hasattr(extract_pairs, 'extract_pairs'))
        except ImportError as e:
            self.fail(f"Could not import pise.extract_pairs: {e}")

    def test_hamming_distance(self):
        from pise.qc import hamming_distance
        self.assertEqual(hamming_distance("CTCTCTAT", "CTCTCTAT"), 0)
        self.assertEqual(hamming_distance("CTCTCTAT", "CTCTCTAA"), 1)
        self.assertEqual(hamming_distance("CTCTCTAT", "TATCCTCT"), 5)

    def test_load_ids_txt_vs_fastq(self):
        from pise.extract_pairs import load_ids
        
        # Test 1: Plain text file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f_txt:
            f_txt.write("id1\n  id2  some_other_info\nid3\n")
            f_txt_path = f_txt.name
            
        try:
            ids = load_ids(f_txt_path)
            self.assertEqual(ids, {"id1", "id2", "id3"})
        finally:
            os.remove(f_txt_path)

        # Test 2: Gzipped FASTQ file
        with tempfile.NamedTemporaryFile(suffix='.fastq.gz', delete=False) as f_fq:
            f_fq_path = f_fq.name
            
        try:
            with gzip.open(f_fq_path, 'wt') as gz_out:
                gz_out.write("@read1 1:N:0:1\nACGT\n+\nIIII\n")
                gz_out.write("@read2\nTGCA\n+\n####\n")
            ids = load_ids(f_fq_path)
            self.assertEqual(ids, {"read1", "read2"})
        finally:
            os.remove(f_fq_path)

    def test_run_preprocess_import(self):
        try:
            from pise.preprocess import run_preprocess
            self.assertTrue(callable(run_preprocess))
        except ImportError as e:
            self.fail(f"Could not import run_preprocess from pise.preprocess: {e}")

if __name__ == '__main__':
    unittest.main()
