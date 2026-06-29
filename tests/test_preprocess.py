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

    def test_extract_pairs_logic(self):
        from pise.extract_pairs import extract_pairs
        
        # Create temp files
        with tempfile.NamedTemporaryFile(suffix='.fastq.gz', delete=False) as f_fw, \
             tempfile.NamedTemporaryFile(suffix='.fastq.gz', delete=False) as f_rv_raw, \
             tempfile.NamedTemporaryFile(suffix='.fastq.gz', delete=False) as f_rv_out:
            
            f_fw_path = f_fw.name
            f_rv_raw_path = f_rv_raw.name
            f_rv_out_path = f_rv_out.name
            
        try:
            # Write 2 matching forward reads
            with gzip.open(f_fw_path, 'wt') as gz_out:
                gz_out.write("@read1 1:N:0:1\nACGT\n+\nIIII\n")
                gz_out.write("@read3\nTGCA\n+\n####\n")
                
            # Write 3 raw reverse reads (read1, read2, read3)
            with gzip.open(f_rv_raw_path, 'wt') as gz_out:
                gz_out.write("@read1 2:N:0:1\nACGT_REV\n+\nIIII_REV\n")
                gz_out.write("@read2\nTGCA_REV\n+\n####_REV\n")
                gz_out.write("@read3\nGGGG\n+\nJJJJ\n")
                
            # Run pairing
            extract_pairs(f_fw_path, f_rv_raw_path, f_rv_out_path)
            
            # Read output
            with gzip.open(f_rv_out_path, 'rt') as gz_in:
                lines = gz_in.readlines()
                
            # Should have read1 and read3, but NOT read2
            self.assertEqual(len(lines), 8) # 2 reads * 4 lines
            self.assertTrue(lines[0].startswith("@read1"))
            self.assertEqual(lines[4], "@read3\n")
        finally:
            os.remove(f_fw_path)
            os.remove(f_rv_raw_path)
            os.remove(f_rv_out_path)

if __name__ == '__main__':
    unittest.main()
