import pytest
import subprocess
import os

# Paths for test files
mock_pairsam_path_dedup = "/workspaces/pairtools/tests/data/mock.4dedup.pairsam"
mock_pairsam_path_dedup_diff_colnames = "/workspaces/pairtools/tests/data/mock.4dedup_diffcolnames.pairsam"

# Temporary output paths
dedup_path = "/tmp/dedup.pairsam"
dups_path = "/tmp/dups.pairsam"
unmapped_path = "/tmp/unmapped.pairsam"

dedup_path_cython = "/tmp/dedup_cython.pairsam"
dups_path_cython = "/tmp/dups_cython.pairsam"
unmapped_path_cython = "/tmp/unmapped_cython.pairsam"

dedup_max_path = "/tmp/dedup_max.pairsam"
dups_max_path = "/tmp/dups_max.pairsam"
unmapped_max_path = "/tmp/unmapped_max.pairsam"

dedup_markdups_path = "/tmp/dedup_markdups.pairsam"
dups_markdups_path = "/tmp/dups_markdups.pairsam"
unmapped_markdups_path = "/tmp/unmapped_markdups.pairsam"

dedup_path_diff_colnames = "/tmp/dedup.diff_colnames.pairsam"
dups_path_diff_colnames = "/tmp/dups.diff_colnames.pairsam"
unmapped_path_diff_colnames = "/tmp/unmapped.diff_colnames.pairsam"

max_mismatch = 1

@pytest.fixture
def setup_dedup():
    try:
        result = subprocess.run(
            [
                "python",
                "-m",
                "pairtools",
                "dedup",
                mock_pairsam_path_dedup,
                "--output",
                dedup_path,
                "--output-dups",
                dups_path,
                "--output-unmapped",
                unmapped_path,
                "--max-mismatch",
                str(max_mismatch),
            ],
            capture_output=True,
            text=True,
        )
        print("Dedup 1 stdout:", result.stdout)
        print("Dedup 1 stderr:", result.stderr)
        result.check_returncode()

        result = subprocess.run(
            [
                "python",
                "-m",
                "pairtools",
                "dedup",
                mock_pairsam_path_dedup,
                "--output",
                dedup_path_cython,
                "--output-dups",
                dups_path_cython,
                "--output-unmapped",
                unmapped_path_cython,
                "--max-mismatch",
                str(max_mismatch),
                "--backend",
                "cython",
            ],
            capture_output=True,
            text=True,
        )
        print("Dedup 2 stdout:", result.stdout)
        print("Dedup 2 stderr:", result.stderr)
        result.check_returncode()

        result = subprocess.run(
            [
                "python",
                "-m",
                "pairtools",
                "dedup",
                mock_pairsam_path_dedup,
                "--output",
                dedup_max_path,
                "--output-dups",
                dups_max_path,
                "--output-unmapped",
                unmapped_max_path,
                "--max-mismatch",
                str(max_mismatch),
                "--method",
                "max",
            ],
            capture_output=True,
            text=True,
        )
        print("Dedup 3 stdout:", result.stdout)
        print("Dedup 3 stderr:", result.stderr)
        result.check_returncode()

        result = subprocess.run(
            [
                "python",
                "-m",
                "pairtools",
                "dedup",
                mock_pairsam_path_dedup,
                "--mark-dups",
                "--output",
                dedup_markdups_path,
                "--output-dups",
                dups_markdups_path,
                "--output-unmapped",
                unmapped_markdups_path,
                "--max-mismatch",
                str(max_mismatch),
            ],
            capture_output=True,
            text=True,
        )
        print("Dedup 4 stdout:", result.stdout)
        print("Dedup 4 stderr:", result.stderr)
        result.check_returncode()

        result = subprocess.run(
            [
                "python",
                "-m",
                "pairtools",
                "dedup",
                mock_pairsam_path_dedup_diff_colnames,
                "--mark-dups",
                "--output",
                dedup_path_diff_colnames,
                "--output-dups",
                dups_path_diff_colnames,
                "--output-unmapped",
                unmapped_path_diff_colnames,
                "--max-mismatch",
                str(max_mismatch),
                "--c1",
                "chrom1",
                "--c2",
                "chrom2",
                "--p1",
                "p1",
                "--p2",
                "p2",
                "--s1",
                "str1",
                "--s2",
                "str2",
            ],
            capture_output=True,
            text=True,
        )
        print("Dedup 5 stdout:", result.stdout)
        print("Dedup 5 stderr:", result.stderr)
        result.check_returncode()

    except subprocess.CalledProcessError as e:
        print(f"Subprocess failed with return code {e.returncode}")
        print("Output:", e.output)
        print("Stderr:", e.stderr)
        raise

def test_mock_pairsam(setup_dedup):
    pairsam_pairs = [
        l.strip().split("\t")
        for l in open(mock_pairsam_path_dedup, "r")
        if not l.startswith("#") and l.strip()
    ]
    for (ddp, up, dp) in [
        (dedup_path, unmapped_path, dups_path),
        (dedup_max_path, unmapped_max_path, dups_max_path),
        (dedup_markdups_path, unmapped_markdups_path, dups_markdups_path),
        (
            dedup_path_diff_colnames,
            unmapped_path_diff_colnames,
            dups_path_diff_colnames,
        ),
    ]:
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
        assert len(dedup_pairs) > 0
        assert len(dup_pairs) > 0