import unittest
import os
import sys

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

if __name__ == '__main__':
    unittest.main()

