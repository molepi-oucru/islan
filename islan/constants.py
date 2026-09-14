# Centralized Constants for the ISLAN Pipeline

# Preprocessing and filtering constants
INDEX_LEN = 8
BATCH_SIZE = 5000
DEFAULT_FULL_OVERLAP_FRACTION = 0.9

# Alignment viewer constants
MAX_ALIGNMENT_READS = 100
GC_WINDOW_SIZE = 21
FLANK_PADDING = 15
PLOT_HEIGHT = 800
READ_PADDING = 6
ZOOM_THRESHOLD = 150
EXTENSION_PADDING = 2500

# Colors for visualization
BASE_COLORS = {
    'A': '#0f9d58', # Green
    'C': '#4285f4', # Blue
    'G': '#f4b400', # Orange/Yellow
    'T': '#db4437'  # Red
}

HEAD_READ_REVERSE_COLOR = 'rgba(147, 197, 253, 0.85)'
HEAD_READ_FORWARD_COLOR = 'rgba(37, 99, 235, 0.85)'
TAIL_READ_REVERSE_COLOR = 'rgba(252, 165, 165, 0.85)'
TAIL_READ_FORWARD_COLOR = 'rgba(220, 38, 38, 0.85)'

HEAD_COV_COLOR = 'rgba(66, 133, 244, 0.8)'
TAIL_COV_COLOR = 'rgba(219, 68, 85, 0.8)'
GC_LINE_COLOR = 'rgba(15, 157, 88, 0.85)'

GENE_FORWARD_COLOR = 'rgba(244, 180, 0, 0.8)'
GENE_REVERSE_COLOR = 'rgba(171, 71, 188, 0.8)'

# Pipeline default parameters
DEFAULT_OUTPUT_DIR = "results_ismapper"
MAPPING_LOG_FILE = "islan_is_mapping.log"
DEFAULT_MIN_CLIP = 10
DEFAULT_MAX_CLIP = 30
DEFAULT_CUTOFF = 6
DEFAULT_MERGING = 100
MAX_PAIRING_DISTANCE = 100
MAX_TSD_OVERLAP = 20
DEFAULT_IS_LENGTH = 4000
DEFAULT_MIN_MAPQ = 30
DEFAULT_FLANK_LEN = 300
DEFAULT_THREADS = 1

# Exhaustive ordered list of detection classes for the HTML summary table.
# Tandem classes are intentionally omitted (rare, not surfaced to end users).
ALL_REPORT_CALL_CLASSES = [
    'known',
    'novel',
    'novel (TSD)',
    'novel (TSD)*',            # negative gap > MAX_TSD_OVERLAP (possible false positive)
    'HEAD-only',
    'TAIL-only',
    'Full Flank Overlap',
]


DEFAULT_MIN_LEN = 20
DEFAULT_I5_MISMATCH = 2
DEFAULT_QC = True
DEFAULT_THREADS_PREPROCESS = 8


class ISElementRegistry:
    """
    Single source of truth for IS element names.
    Parses targets.fasta (or primers.fasta) dynamically — no hardcoded names.

    Name formats derived from FASTA header '>ST16:ISKpn26_IS5:HEAD':
      full_name:  'ISKpn26_IS5'  -> used by load_primers(), load_known_indices()
      short_name: 'ISKpn26'      -> used for demux categories, QC charts, --is_name

    Usage::

        registry = ISElementRegistry('/path/to/targets.fasta')
        full, short = registry.resolve('502-ISKpn26')
        # -> ('ISKpn26_IS5', 'ISKpn26')
    """

    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.short_to_full = {}  # e.g. {'ISKpn26': 'ISKpn26_IS5', ...}
        self.full_to_short = {}  # e.g. {'ISKpn26_IS5': 'ISKpn26', ...}
        self._parse(fasta_path)

    def _parse(self, fasta_path):
        """Parse a FASTA file and extract unique IS element names."""
        from Bio import SeqIO
        seen_full = set()
        try:
            with open(fasta_path) as f:
                for record in SeqIO.parse(f, 'fasta'):
                    parts = record.id.split(':')
                    if len(parts) < 2:
                        continue
                    full = parts[1]             # e.g. 'ISKpn26_IS5'
                    short = full.split('_')[0]  # e.g. 'ISKpn26'
                    if full not in seen_full:
                        seen_full.add(full)
                        self.short_to_full[short] = full
                        self.full_to_short[full] = short
        except (FileNotFoundError, IOError):
            # Graceful degradation: registry stays empty, callers should warn
            pass

    @property
    def all_full_names(self):
        """Ordered list of full IS names, e.g. ['IS1R_IS1', 'ISKpn26_IS5']."""
        return list(self.full_to_short.keys())

    @property
    def all_short_names(self):
        """Ordered list of short IS names, e.g. ['IS1R', 'ISKpn26']."""
        return list(self.short_to_full.keys())

    def resolve(self, name):
        """
        Fuzzy-match a sample/filename string to (full_name, short_name).

        Matches by checking whether any known short name is a case-insensitive
        substring of the input, e.g.::

            '502-ISKpn26' -> ('ISKpn26_IS5', 'ISKpn26')
            '278-IS1R'    -> ('IS1R_IS1',    'IS1R')

        Returns (None, None) if no match is found.
        """
        name_upper = name.upper()
        for short, full in self.short_to_full.items():
            if short.upper() in name_upper:
                return full, short
        return None, None
