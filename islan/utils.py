import sys
import logging
from Bio import SeqIO
from Bio.Seq import Seq

def load_primers(fasta_file_path, target_is_element):
    """
    Load primers for a specific IS element from primers.fasta.
    Returns a dictionary containing 'HEAD', 'TAIL', 'P_UP', and 'P_DOWN'.
    Note: P_UP and P_DOWN trimmed of their 5' 8bp index sequence.
    """
    primers = {}
    try:
        with open(fasta_file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                parts = record.id.split(':')
                if len(parts) >= 2 and parts[1] == target_is_element:
                    if ":HEAD" in record.id:
                        primers['HEAD'] = str(Seq(record.seq).reverse_complement())
                    elif ":TAIL" in record.id:
                        primers['TAIL'] = str(record.seq)
                    elif ":P_UP" in record.id:
                        primers['P_UP'] = str(record.seq)[8:]
                    elif ":P_DOWN" in record.id:
                        primers['P_DOWN'] = str(record.seq)[8:]
    except FileNotFoundError:
        logging.error(f"Primer FASTA file not found at {fasta_file_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An error occurred while loading primers: {e}")
        sys.exit(1)
    return primers

def load_known_indices(primers_path):
    """
    Load known 8bp i5 indices for each IS element from primers.fasta.
    """
    known_indices = {}
    try:
        with open(primers_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id.endswith(":P_UP"):
                    parts = record.id.split(':')
                    if len(parts) >= 2:
                        is_element = parts[1]
                        known_indices[is_element] = str(record.seq[:8])
    except Exception as e:
        logging.warning(f"Could not load known indices from {primers_path}: {e}")
    return known_indices
