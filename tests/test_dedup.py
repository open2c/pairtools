import pytest
import subprocess
import os

# Dynamically determine the test directory
testdir = os.path.dirname(os.path.realpath(__file__))

# Paths for test input files (relative to testdir)
mock_pairsam_path_dedup = os.path.join(testdir, "data", "mock.4dedup.pairsam")
mock_pairsam_path_dedup_diff_colnames = os.path.join(testdir, "data", "mock.4dedup_diffcolnames.pairsam")

# Global variables for output paths
dedup_path = None
dups_path = None
unmapped_path = None
dedup_path_cython = None
dups_path_cython = None
unmapped_path_cython = None
dedup_max_path = None
dups_max_path = None
unmapped_max_path = None
dedup_markdups_path = None
dups_markdups_path = None
unmapped_markdups_path = None
dedup_path_diff_colnames = None
dups_path_diff_colnames = None
unmapped_path_diff_colnames = None

# Constant for max mismatch
max_mismatch = 1

@pytest.fixture
def setup_dedup(tmp_path):
    """Fixture to set up deduplication by running pairtools dedup commands."""
    # Declare global variables to modify them within the fixture
    global dedup_path, dups_path, unmapped_path, dedup_path_cython, dups_path_cython, unmapped_path_cython
    global dedup_max_path, dups_max_path, unmapped_max_path, dedup_markdups_path, dups_markdups_path, unmapped_markdups_path
    global dedup_path_diff_colnames, dups_path_diff_colnames, unmapped_path_diff_colnames

    # Define output paths within tmp_path
    dedup_path = str(tmp_path / "dedup.pairsam")
    dups_path = str(tmp_path / "dups.pairsam")
    unmapped_path = str(tmp_path / "unmapped.pairsam")
    dedup_path_cython = str(tmp_path / "dedup_cython.pairsam")
    dups_path_cython = str(tmp_path / "dups_cython.pairsam")
    unmapped_path_cython = str(tmp_path / "unmapped_cython.pairsam")
    dedup_max_path = str(tmp_path / "dedup_max.pairsam")
    dups_max_path = str(tmp_path / "dups_max.pairsam")
    unmapped_max_path = str(tmp_path / "unmapped_max.pairsam")
    dedup_markdups_path = str(tmp_path / "dedup_markdups.pairsam")
    dups_markdups_path = str(tmp_path / "dups_markdups.pairsam")
    unmapped_markdups_path = str(tmp_path / "unmapped_markdups.pairsam")
    dedup_path_diff_colnames = str(tmp_path / "dedup.diff_colnames.pairsam")
    dups_path_diff_colnames = str(tmp_path / "dups.diff_colnames.pairsam")
    unmapped_path_diff_colnames = str(tmp_path / "unmapped.diff_colnames.pairsam")

    # List of pairtools dedup commands to execute
    commands = [
        # Standard dedup
        [
            "python", "-m", "pairtools", "dedup",
            mock_pairsam_path_dedup,
            "--output", dedup_path,
            "--output-dups", dups_path,
            "--output-unmapped", unmapped_path,
            "--max-mismatch", str(max_mismatch),
            "--c1", "chrom1", "--p1", "pos1", "--s1", "strand1",
            "--c2", "chrom2", "--p2", "pos2", "--s2", "strand2",
        ],
        # Dedup with Cython backend
        [
            "python", "-m", "pairtools", "dedup",
            mock_pairsam_path_dedup,
            "--output", dedup_path_cython,
            "--output-dups", dups_path_cython,
            "--output-unmapped", unmapped_path_cython,
            "--max-mismatch", str(max_mismatch),
            "--backend", "cython",
            "--c1", "chrom1", "--p1", "pos1", "--s1", "strand1",
            "--c2", "chrom2", "--p2", "pos2", "--s2", "strand2",
        ],
        # Dedup with max method
        [
            "python", "-m", "pairtools", "dedup",
            mock_pairsam_path_dedup,
            "--output", dedup_max_path,
            "--output-dups", dups_max_path,
            "--output-unmapped", unmapped_max_path,
            "--max-mismatch", str(max_mismatch),
            "--method", "max",
            "--c1", "chrom1", "--p1", "pos1", "--s1", "strand1",
            "--c2", "chrom2", "--p2", "pos2", "--s2", "strand2",
        ],
        # Dedup with mark-dups
        [
            "python", "-m", "pairtools", "dedup",
            mock_pairsam_path_dedup,
            "--mark-dups",
            "--output", dedup_markdups_path,
            "--output-dups", dups_markdups_path,
            "--output-unmapped", unmapped_markdups_path,
            "--max-mismatch", str(max_mismatch),
            "--c1", "chrom1", "--p1", "pos1", "--s1", "strand1",
            "--c2", "chrom2", "--p2", "pos2", "--s2", "strand2",
        ],
        # Dedup with different column names
        [
            "python", "-m", "pairtools", "dedup",
            mock_pairsam_path_dedup_diff_colnames,
            "--mark-dups",
            "--output", dedup_path_diff_colnames,
            "--output-dups", dups_path_diff_colnames,
            "--output-unmapped", unmapped_path_diff_colnames,
            "--max-mismatch", str(max_mismatch),
            "--c1", "chrom1", "--c2", "chrom2",
            "--p1", "p1", "--p2", "p2",
            "--s1", "str1", "--s2", "str2",
        ],
    ]

    # Execute each command
    for i, cmd in enumerate(commands, 1):
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )
        print(f"Dedup {i} stdout:", result.stdout)
        print(f"Dedup {i} stderr:", result.stderr)
        result.check_returncode()  # Raise an exception if the command fails

def test_mock_pairsam(setup_dedup):
    """Test that deduplicated output files contain data."""
    # Read the input pairsam file
    pairsam_pairs = [
        l.strip().split("\t")
        for l in open(mock_pairsam_path_dedup, "r")
        if not l.startswith("#") and l.strip()
    ]

    # List of output path tuples to check
    output_paths = [
        (dedup_path, unmapped_path, dups_path),
        (dedup_path_cython, unmapped_path_cython, dups_path_cython),
        (dedup_max_path, unmapped_max_path, dups_max_path),
        (dedup_markdups_path, unmapped_markdups_path, dups_markdups_path),
        (dedup_path_diff_colnames, unmapped_path_diff_colnames, dups_path_diff_colnames),
    ]

    # Check each set of output files
    for ddp, up, dp in output_paths:
        dedup_pairs = [
            l.strip().split("\t")
            for l in open(ddp, "r")
            if not l.startswith("#") and l.strip()
        ]
        unmapped_pairs = [
            l.strip().split("\t")
            for l in open(up, "r")
            if not l.startswith("#") and l.strip()
        ]
        dup_pairs = [
            l.strip().split("\t")
            for l in open(dp, "r")
            if not l.startswith("#") and l.strip()
        ]
        assert len(dedup_pairs) > 0, f"No deduplicated pairs found in {ddp}"
        assert len(dup_pairs) > 0, f"No duplicate pairs found in {dp}"