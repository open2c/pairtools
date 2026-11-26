# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import pytest
from pairtools.lib import pairsam_format

testdir = os.path.dirname(os.path.realpath(__file__))
mock_pairsam_path = os.path.join(testdir, "data", "mock.pairsam")
mock_chromsizes_path = os.path.join(testdir, "data", "mock.chrom.sizes")


def test_preserve():
    try:
        result = subprocess.check_output(
            ["python", "-m", "pairtools", "select", "True", mock_pairsam_path],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    pairsam_body = [
        l.strip()
        for l in open(mock_pairsam_path, "r")
        if not l.startswith("#") and l.strip()
    ]
    output_body = [
        l.strip() for l in result.split("\n") if not l.startswith("#") and l.strip()
    ]
    assert all(l in pairsam_body for l in output_body)


def test_equal():
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "select",
                '(pair_type == "RU") or (pair_type == "UR") or (pair_type == "UU")',
                mock_pairsam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e
    print(result)

    pairsam_body = [
        l.strip()
        for l in open(mock_pairsam_path, "r")
        if not l.startswith("#") and l.strip()
    ]
    output_body = [
        l.strip() for l in result.split("\n") if not l.startswith("#") and l.strip()
    ]

    assert all(l.split("\t")[7] in ["RU", "UR", "UU"] for l in output_body)
    assert all(
        l in output_body for l in pairsam_body if l.split("\t")[7] in ["RU", "UR", "UU"]
    )


def test_csv():
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "select",
                'csv_match(pair_type, "RU,UR,UU")',
                mock_pairsam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e
    print(result)

    pairsam_body = [
        l.strip()
        for l in open(mock_pairsam_path, "r")
        if not l.startswith("#") and l.strip()
    ]
    output_body = [
        l.strip() for l in result.split("\n") if not l.startswith("#") and l.strip()
    ]

    assert all(l.split("\t")[7] in ["RU", "UR", "UU"] for l in output_body)
    assert all(
        l in output_body for l in pairsam_body if l.split("\t")[7] in ["RU", "UR", "UU"]
    )


def test_wildcard():
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "select",
                'wildcard_match(pair_type, "*U")',
                mock_pairsam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e
    print(result)

    pairsam_body = [
        l.strip()
        for l in open(mock_pairsam_path, "r")
        if not l.startswith("#") and l.strip()
    ]
    output_body = [
        l.strip() for l in result.split("\n") if not l.startswith("#") and l.strip()
    ]

    assert all(l.split("\t")[7] in ["NU", "MU", "RU", "UU"] for l in output_body)
    assert all(
        l in output_body
        for l in pairsam_body
        if l.split("\t")[7] in ["NU", "MU", "RU", "UU"]
    )


def test_regex():
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "select",
                'regex_match(pair_type, "[NM]U")',
                mock_pairsam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e
    print(result)

    pairsam_body = [
        l.strip()
        for l in open(mock_pairsam_path, "r")
        if not l.startswith("#") and l.strip()
    ]
    output_body = [
        l.strip() for l in result.split("\n") if not l.startswith("#") and l.strip()
    ]

    assert all(l.split("\t")[7] in ["NU", "MU"] for l in output_body)
    assert all(
        l in output_body for l in pairsam_body if l.split("\t")[7] in ["NU", "MU"]
    )


def test_chrom_subset():
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "select",
                "True",
                "--chrom-subset",
                mock_chromsizes_path,
                mock_pairsam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    pairsam_body = [
        l.strip()
        for l in open(mock_pairsam_path, "r")
        if not l.startswith("#") and l.strip()
    ]
    output_body = [
        l.strip() for l in result.split("\n") if not l.startswith("#") and l.strip()
    ]
    output_header = [
        l.strip() for l in result.split("\n") if l.startswith("#") and l.strip()
    ]

    chroms_from_chrom_field = [
        l.strip().split()[1:]
        for l in result.split("\n")
        if l.startswith("#chromosomes:")
    ][0]

    assert set(chroms_from_chrom_field) == set(["chr1", "chr2"])

    chroms_from_chrom_sizes = [
        l.strip().split()[1] for l in result.split("\n") if l.startswith("#chromsize:")
    ]

    assert set(chroms_from_chrom_sizes) == set(["chr1", "chr2"])


def test_remove_columns():
    """Test removal of columns from the file
    Example run:
    pairtools select True --remove-columns sam1,sam2 tests/data/mock.pairsam
    """

    mock_pairs_path = os.path.join(testdir, "data", "mock.pairsam")
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "select",
                "True",
                "--remove-columns",
                "sam1,sam2",
                mock_pairs_path,
            ],
        ).decode("ascii")

    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e

    # check if the columns are removed properly:
    pairsam_header = [l.strip() for l in result.split("\n") if l.startswith("#")]
    for l in pairsam_header:
        if l.startswith("#columns:"):
            line = l.strip()
            assert (
                line
                == "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"
            )

    # check that the pairs got assigned properly
    for l in result.split("\n"):
        if l.startswith("#") or not l:
            continue

        assert len(l.split(pairsam_format.PAIRSAM_SEP)) == 8


def test_region_match():
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "select",
                'region_match(chrom1, pos1, "chr1", 0, 50)',
                mock_pairsam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e
    print(result)

    pairsam_body = [
        l.strip()
        for l in open(mock_pairsam_path, "r")
        if not l.startswith("#") and l.strip()
    ]
    output_body = [
        l.strip() for l in result.split("\n") if not l.startswith("#") and l.strip()
    ]

    # Verify all output rows have chrom1="chr1" and pos1 within range
    for l in output_body:
        fields = l.split("\t")
        chrom1, pos1 = fields[1], int(fields[2])
        assert chrom1 == "chr1"
        assert 0 <= pos1 <= 50

    # Verify all matching rows from input are in output
    for l in pairsam_body:
        fields = l.split("\t")
        chrom1, pos1 = fields[1], int(fields[2])
        if chrom1 == "chr1" and 0 <= pos1 <= 50:
            assert l in output_body


def test_region_match_no_end():
    try:
        result = subprocess.check_output(
            [
                "python",
                "-m",
                "pairtools",
                "select",
                'region_match(chrom1, pos1, "chr1", 50)',
                mock_pairsam_path,
            ],
        ).decode("ascii")
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(sys.exc_info())
        raise e
    print(result)

    pairsam_body = [
        l.strip()
        for l in open(mock_pairsam_path, "r")
        if not l.startswith("#") and l.strip()
    ]
    output_body = [
        l.strip() for l in result.split("\n") if not l.startswith("#") and l.strip()
    ]

    # Verify all output rows have chrom1="chr1" and pos1 >= 50
    for l in output_body:
        fields = l.split("	")
        chrom1, pos1 = fields[1], int(fields[2])
        assert chrom1 == "chr1"
        assert pos1 >= 50

    # Verify all matching rows from input are in output
    for l in pairsam_body:
        fields = l.split("	")
        chrom1, pos1 = fields[1], int(fields[2])
        if chrom1 == "chr1" and pos1 >= 100:
            assert l in output_body


EXTRA_COLS_PAIRS = (
    "## pairs format v1.0.0\n"
    "#shape: upper triangle\n"
    "#genome_assembly: unknown\n"
    "#chromosomes: chr1 chr2\n"
    "#chromsize: chr1 1000\n"
    "#chromsize: chr2 1000\n"
    "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type "
    "walk_pair_index walk_pair_type read_len1 read_len2 mapq1 mapq2\n"
    "r1\tchr1\t10\tchr1\t100\t+\t+\tUU\t0\tR1\t100\t100\t60\t60\n"
    "r2\tchr1\t20\tchr1\t200\t+\t+\tUU\t1\tR2\t100\t100\t0\t0\n"
    "r3\tchr1\t30\tchr2\t300\t+\t-\tUU\t0\tR1\t100\t100\t60\t5\n"
    "r4\tchr2\t40\tchr2\t400\t-\t+\tUU\t2\tR2\t50\t120\t45\t55\n"
)


def _parse_pairs_output(text):
    header = [line for line in text.splitlines() if line.startswith("#")]
    body = [line for line in text.splitlines() if line and not line.startswith("#")]
    cols_lines = [line for line in header if line.startswith("#columns:")]
    assert len(cols_lines) == 1
    columns = cols_lines[0].split()[1:]
    rows = [line.split(pairsam_format.PAIRSAM_SEP) for line in body]
    return columns, rows


def test_remove_columns_filter_after_removed_and_output_rest(tmp_path):
    """Regression test for #279: filtering on columns after removed columns."""
    inp = tmp_path / "extra.pairs"
    inp.write_text(EXTRA_COLS_PAIRS)
    out_rest = tmp_path / "rest.pairs"

    condition = "(mapq1>=30) and (mapq2>=30)"

    ctrl_out = subprocess.check_output(
        [sys.executable, "-m", "pairtools", "select", condition, str(inp)],
    ).decode("ascii")
    _, ctrl_rows = _parse_pairs_output(ctrl_out)
    ctrl_read_ids = [r[0] for r in ctrl_rows]
    assert ctrl_read_ids == ["r1", "r4"]

    out = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "read_len1,read_len2",
            "--output-rest",
            str(out_rest),
            condition,
            str(inp),
        ],
    ).decode("ascii")

    expected_cols = [
        "readID",
        "chrom1",
        "pos1",
        "chrom2",
        "pos2",
        "strand1",
        "strand2",
        "pair_type",
        "walk_pair_index",
        "walk_pair_type",
        "mapq1",
        "mapq2",
    ]

    sel_cols, sel_rows = _parse_pairs_output(out)
    assert sel_cols == expected_cols
    assert [r[0] for r in sel_rows] == ctrl_read_ids
    assert all(len(r) == len(expected_cols) for r in sel_rows)
    assert [(r[-2], r[-1]) for r in sel_rows] == [("60", "60"), ("45", "55")]

    rest_cols, rest_rows = _parse_pairs_output(out_rest.read_text())
    assert rest_cols == expected_cols
    assert [r[0] for r in rest_rows] == ["r2", "r3"]
    assert all(len(r) == len(expected_cols) for r in rest_rows)
    assert [(r[-2], r[-1]) for r in rest_rows] == [("0", "0"), ("60", "5")]


def test_remove_columns_filter_on_removed_column(tmp_path):
    """Filtering on a column that is itself removed should evaluate against input rows."""
    inp = tmp_path / "extra.pairs"
    inp.write_text(EXTRA_COLS_PAIRS)

    out = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "read_len1",
            "-t",
            "read_len1",
            "int",
            "(read_len1>=80) and (mapq1>=30)",
            str(inp),
        ],
    ).decode("ascii")

    cols, rows = _parse_pairs_output(out)
    assert "read_len1" not in cols
    assert len(cols) == 13
    assert [r[0] for r in rows] == ["r1", "r3"]
    assert all(len(r) == 13 for r in rows)


def test_remove_columns_with_type_cast_kept_and_removed(tmp_path):
    """--type-cast must index original input columns for both kept and removed columns."""
    pairs_text = (
        "## pairs format v1.0.0\n"
        "#shape: upper triangle\n"
        "#genome_assembly: unknown\n"
        "#chromsize: chr1 1000\n"
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type "
        "read_len1 read_len2 score mapq1 mapq2\n"
        "r1\tchr1\t10\tchr1\t100\t+\t+\tUU\t100\t100\t0.9\t60\t60\n"
        "r2\tchr1\t20\tchr1\t200\t+\t+\tUU\t100\t100\t0.2\t60\t60\n"
        "r3\tchr1\t30\tchr1\t300\t+\t+\tUU\t50\t100\t0.8\t60\t60\n"
        "r4\tchr1\t40\tchr1\t400\t+\t+\tUU\t100\t100\t0.9\t10\t10\n"
    )
    inp = tmp_path / "scores.pairs"
    inp.write_text(pairs_text)

    # Type-cast on a kept column appearing after removed columns
    out_kept = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "read_len1,read_len2",
            "-t",
            "score",
            "float",
            "(score>=0.5) and (mapq1>=30)",
            str(inp),
        ],
    ).decode("ascii")
    cols_kept, rows_kept = _parse_pairs_output(out_kept)
    assert "read_len1" not in cols_kept and "read_len2" not in cols_kept
    assert "score" in cols_kept
    assert [r[0] for r in rows_kept] == ["r1", "r3"]
    assert all(len(r) == len(cols_kept) for r in rows_kept)

    # Type-cast on both a removed column and a kept column
    out_both = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "read_len1,score",
            "-t",
            "read_len1",
            "int",
            "-t",
            "score",
            "float",
            "(read_len1>=80) and (score>=0.5) and (mapq1>=30)",
            str(inp),
        ],
    ).decode("ascii")
    cols_both, rows_both = _parse_pairs_output(out_both)
    assert "read_len1" not in cols_both and "score" not in cols_both
    assert [r[0] for r in rows_both] == ["r1"]
    assert all(len(r) == len(cols_both) for r in rows_both)


def test_remove_columns_non_contiguous(tmp_path):
    """Removing non-contiguous extra columns while filtering on interleaved kept columns."""
    inp = tmp_path / "extra.pairs"
    inp.write_text(EXTRA_COLS_PAIRS)

    # Extra cols: walk_pair_index (1st), walk_pair_type (2nd), read_len1 (3rd), read_len2 (4th)
    # Remove 1st and 3rd; filter on 2nd, 4th, and mapq2
    condition = '(walk_pair_type == "R2") and (read_len2 >= 110) and (mapq2 >= 30)'

    ctrl_out = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "-t",
            "read_len2",
            "int",
            condition,
            str(inp),
        ],
    ).decode("ascii")
    _, ctrl_rows = _parse_pairs_output(ctrl_out)
    assert [r[0] for r in ctrl_rows] == ["r4"]

    out = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "walk_pair_index,read_len1",
            "-t",
            "read_len2",
            "int",
            condition,
            str(inp),
        ],
    ).decode("ascii")
    cols, rows = _parse_pairs_output(out)
    assert cols == [
        "readID",
        "chrom1",
        "pos1",
        "chrom2",
        "pos2",
        "strand1",
        "strand2",
        "pair_type",
        "walk_pair_type",
        "read_len2",
        "mapq1",
        "mapq2",
    ]
    assert [r[0] for r in rows] == [r[0] for r in ctrl_rows]
    assert rows[0][8:] == ["R2", "120", "45", "55"]


def test_remove_columns_duplicate_and_chrom_subset(tmp_path):
    """Duplicate names in --remove-columns and --remove-columns with --chrom-subset."""
    inp = tmp_path / "extra.pairs"
    inp.write_text(EXTRA_COLS_PAIRS)
    chroms = tmp_path / "chroms.txt"
    chroms.write_text("chr1\t1000\n")

    # Duplicate column name in --remove-columns
    out_dup = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "read_len1,read_len1",
            "(mapq1>=30) and (mapq2>=30)",
            str(inp),
        ],
    ).decode("ascii")
    cols_dup, rows_dup = _parse_pairs_output(out_dup)
    assert "read_len1" not in cols_dup
    assert [r[0] for r in rows_dup] == ["r1", "r4"]
    assert all(len(r) == len(cols_dup) for r in rows_dup)

    # Removing a core column before chrom1/chrom2 alongside --chrom-subset
    out_chrom = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "readID",
            "--chrom-subset",
            str(chroms),
            "(mapq1>=30) and (mapq2>=30)",
            str(inp),
        ],
    ).decode("ascii")
    cols_chrom, rows_chrom = _parse_pairs_output(out_chrom)
    assert cols_chrom[0] == "chrom1"
    assert len(rows_chrom) == 1
    assert rows_chrom[0][0] == "chr1" and rows_chrom[0][1] == "10"
    assert all(len(r) == len(cols_chrom) for r in rows_chrom)

    # Nonexistent column name in --remove-columns has no effect on header or rows
    out_invalid = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "nonexistent_col",
            "(mapq1>=30) and (mapq2>=30)",
            str(inp),
        ],
    ).decode("ascii")
    cols_inv, rows_inv = _parse_pairs_output(out_invalid)
    assert len(cols_inv) == 14
    assert [r[0] for r in rows_inv] == ["r1", "r4"]
    assert all(len(r) == 14 for r in rows_inv)

    # File lacking #columns: header with --remove-columns must not wipe body rows
    no_cols_pairs = (
        "## pairs format v1.0.0\n"
        "#chromsize: chr1 1000\n"
        "r1\tchr1\t10\tchr1\t100\t+\t+\tUU\t.\t.\n"
    )
    no_cols_inp = tmp_path / "no_cols.pairs"
    no_cols_inp.write_text(no_cols_pairs)
    out_nocols = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "sam1",
            "pos1 >= 10",
            str(no_cols_inp),
        ],
    ).decode("ascii")
    body_nocols = [
        line for line in out_nocols.splitlines() if line and not line.startswith("#")
    ]
    assert body_nocols == ["r1\tchr1\t10\tchr1\t100\t+\t+\tUU\t.\t."]


def test_evaluate_stream_does_not_mutate_global_types():
    """Calling evaluate_stream with type_cast must not leak into module-level TYPES."""
    from pairtools.lib import select as select_lib

    original_types = select_lib.TYPES.copy()
    stream = ["r1\tchr1\t10\tchr1\t20\t+\t+\tUU\t100\n"]
    cols = [
        "readID",
        "chrom1",
        "pos1",
        "chrom2",
        "pos2",
        "strand1",
        "strand2",
        "pair_type",
        "custom_num",
    ]
    res = list(
        select_lib.evaluate_stream(
            stream, "custom_num >= 50", cols, type_cast=(("custom_num", "int"),)
        )
    )
    assert res == [(True, stream[0])]
    assert select_lib.TYPES == original_types


def test_select_preserves_quoted_column_names(tmp_path):
    """Quoted string literals matching column names must not be substituted."""
    pairs_text = (
        "## pairs format v1.0.0\n"
        "#shape: upper triangle\n"
        "#genome_assembly: unknown\n"
        "#chromsize: chr1 1000\n"
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type "
        "read_len1 mapq1\n"
        "read_len1\tchr1\t10\tchr1\t100\t+\t+\tUU\t80\t60\n"
        "mapq1\tchr1\t20\tchr1\t200\t+\t+\tUU\t100\t40\n"
        "r3\tchr1\t30\tchr1\t300\t+\t+\tUU\t80\t20\n"
    )
    inp = tmp_path / "literals.pairs"
    inp.write_text(pairs_text)
    out_rest = tmp_path / "rest.pairs"

    # Double-quoted literal matching a removed column name (referee counterexample)
    out = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "read_len1",
            "--output-rest",
            str(out_rest),
            'readID=="read_len1"',
            str(inp),
        ],
    ).decode("ascii")
    cols, rows = _parse_pairs_output(out)
    assert "read_len1" not in cols
    assert [r[0] for r in rows] == ["read_len1"]
    assert all(len(r) == len(cols) for r in rows)

    rest_cols, rest_rows = _parse_pairs_output(out_rest.read_text())
    assert rest_cols == cols
    assert [r[0] for r in rest_rows] == ["mapq1", "r3"]

    # Single-quoted literal matching a removed column name
    out_single = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "read_len1",
            "readID=='read_len1'",
            str(inp),
        ],
    ).decode("ascii")
    _, rows_single = _parse_pairs_output(out_single)
    assert [r[0] for r in rows_single] == ["read_len1"]

    # Without --remove-columns, matching a kept column name ("mapq1")
    out_kept = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            'readID=="mapq1"',
            str(inp),
        ],
    ).decode("ascii")
    _, rows_kept = _parse_pairs_output(out_kept)
    assert [r[0] for r in rows_kept] == ["mapq1"]

    # Combined literal and identifier for the same column name, plus raw/triple/bytes
    out_combo = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "read_len1",
            '(readID=="read_len1") and (read_len1=="80") '
            'and (readID==r"read_len1") and (readID=="""read_len1""") '
            'and (readID.encode()==b"read_len1")',
            str(inp),
        ],
    ).decode("ascii")
    cols_combo, rows_combo = _parse_pairs_output(out_combo)
    assert "read_len1" not in cols_combo
    assert [r[0] for r in rows_combo] == ["read_len1"]


def test_select_preserves_attributes_keyword_args_and_syntax_errors(tmp_path):
    """Attribute accesses and keyword argument names must not be rewritten."""
    from pairtools.lib import select as select_lib

    pairs_text = (
        "## pairs format v1.0.0\n"
        "#shape: upper triangle\n"
        "#genome_assembly: unknown\n"
        "#chromsize: chr1 1000\n"
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type "
        "read_len1 replace\n"
        "read_len1\tchr1\t10\tchr1\t100\t+\t+\tUU\t80\tfoo\n"
        "other\tchr1\t20\tchr1\t200\t+\t+\tUU\t50\tbar\n"
    )
    inp = tmp_path / "attr_kw.pairs"
    inp.write_text(pairs_text)

    # Attribute access (str.replace) + literal ("read_len1") + column identifiers
    out_attr = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            '(readID.replace("read_len1", "ok") == "ok") '
            'and (read_len1 == "80") and (replace == "foo")',
            str(inp),
        ],
    ).decode("ascii")
    _, rows_attr = _parse_pairs_output(out_attr)
    assert [r[0] for r in rows_attr] == ["read_len1"]

    # Keyword argument name colliding with a column name via --startup-code
    out_kw = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--startup-code",
            'def check_len(read_len1, val): return read_len1 == 80 and val == "80"',
            'check_len(read_len1=80, val=read_len1) and (readID == "read_len1")',
            str(inp),
        ],
    ).decode("ascii")
    _, rows_kw = _parse_pairs_output(out_kw)
    assert [r[0] for r in rows_kw] == ["read_len1"]

    # Malformed syntax (TokenError in tokenizer) still raises SyntaxError cleanly
    cols = [
        "readID",
        "chrom1",
        "pos1",
        "chrom2",
        "pos2",
        "strand1",
        "strand2",
        "pair_type",
    ]
    with pytest.raises(SyntaxError):
        list(
            select_lib.evaluate_stream(
                ["r1\tchr1\t10\tchr1\t20\t+\t+\tUU\n"],
                'readID == "unterminated',
                cols,
            )
        )


def test_select_fstrings_preserve_debug_prefix_and_bind_columns(tmp_path):
    """F-strings bind column names cleanly and preserve f'{col=}' debug text."""
    import pandas as pd
    from pairtools.lib import select as select_lib

    pairs_text = (
        "## pairs format v1.0.0\n"
        "#shape: upper triangle\n"
        "#genome_assembly: unknown\n"
        "#chromsize: chr1 1000\n"
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type "
        "read_len1 read_len2 mapq1 mapq2\n"
        "r1\tchr1\t10\tchr1\t100\t+\t+\tUU\t100\t100\t30\t60\n"
        "r2\tchr1\t20\tchr1\t200\t+\t+\tUU\t100\t100\t0\t60\n"
        "r3\tchr1\t30\tchr1\t300\t+\t+\tUU\t30\t100\t60\t5\n"
    )
    inp = tmp_path / "fstrings.pairs"
    inp.write_text(pairs_text)

    # f"{mapq1}" == "30" without --remove-columns
    out_plain = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            'f"{mapq1}" == "30"',
            str(inp),
        ],
    ).decode("ascii")
    _, rows_plain = _parse_pairs_output(out_plain)
    assert [r[0] for r in rows_plain] == ["r1"]

    # f"{mapq1}" == "30" with --remove-columns read_len1,read_len2
    out_rm = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "read_len1,read_len2",
            'f"{mapq1}" == "30"',
            str(inp),
        ],
    ).decode("ascii")
    cols_rm, rows_rm = _parse_pairs_output(out_rm)
    assert "read_len1" not in cols_rm and "read_len2" not in cols_rm
    assert [r[0] for r in rows_rm] == ["r1"]
    assert all(len(r) == len(cols_rm) for r in rows_rm)

    # Column binding in f-string beats a same-named startup variable
    out_startup = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--startup-code",
            "mapq1 = 30",
            'f"{mapq1}" == "30"',
            str(inp),
        ],
    ).decode("ascii")
    _, rows_startup = _parse_pairs_output(out_startup)
    assert [r[0] for r in rows_startup] == ["r1"]

    # f"{mapq1=}" debug prefix must stay "mapq1=", not "int(COLS[...])="
    out_debug = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "--remove-columns",
            "read_len1,read_len2",
            'f"{mapq1=}" == "mapq1=30"',
            str(inp),
        ],
    ).decode("ascii")
    cols_debug, rows_debug = _parse_pairs_output(out_debug)
    assert "read_len1" not in cols_debug and "read_len2" not in cols_debug
    assert [r[0] for r in rows_debug] == ["r1"]

    # Non-pandas evaluate_df also preserves f"{mapq1=}" debug text
    df = pd.DataFrame({"readID": ["r1", "r2"], "mapq1": [30, 0]})
    mask = select_lib.evaluate_df(df, 'f"{mapq1=}" == "mapq1=30"', engine="python")
    assert mask == [True, False]


def test_select_scoping_lambda_comprehensions_and_walrus(tmp_path):
    """Lambda params, comprehension targets, and walrus targets are not clobbered."""
    pairs_text = (
        "## pairs format v1.0.0\n"
        "#shape: upper triangle\n"
        "#genome_assembly: unknown\n"
        "#chromsize: chr1 1000\n"
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type "
        "read_len1 read_len2 mapq1 mapq2\n"
        "r1\tchr1\t10\tchr1\t100\t+\t+\tUU\t100\t100\t30\t60\n"
        "r2\tchr1\t20\tchr1\t200\t+\t+\tUU\t100\t100\t0\t60\n"
        "r3\tchr1\t30\tchr1\t300\t+\t+\tUU\t30\t100\t60\t5\n"
    )
    inp = tmp_path / "scoping.pairs"
    inp.write_text(pairs_text)

    # Shadowed lambda parameter and enclosing-scope lambda default
    out_lambda = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "(lambda mapq1: mapq1 == 1)(1) and (lambda mapq1=mapq1: mapq1 == 30)()",
            str(inp),
        ],
    ).decode("ascii")
    _, rows_lambda = _parse_pairs_output(out_lambda)
    assert [r[0] for r in rows_lambda] == ["r1"]

    # Comprehensions: shadowed target vs outer iterable bound to column
    out_comp = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "([mapq1 for mapq1 in [1, 2]] == [1, 2]) "
            "and ([x for x in [mapq1]] == [30]) "
            "and ([mapq1 for mapq1 in [mapq1]] == [30]) "
            "and ({mapq1 for mapq1 in [1, 2]} == {1, 2}) "
            "and ({mapq1: mapq2 for mapq1, mapq2 in [(1, 2)]} == {1: 2}) "
            "and (list(mapq1 for mapq1 in [1, 2]) == [1, 2])",
            str(inp),
        ],
    ).decode("ascii")
    _, rows_comp = _parse_pairs_output(out_comp)
    assert [r[0] for r in rows_comp] == ["r1"]

    # Walrus expressions: value binds column; target is Store() and not replaced
    out_walrus_val = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "(v := mapq1) >= 30 and v == mapq1",
            str(inp),
        ],
    ).decode("ascii")
    _, rows_walrus_val = _parse_pairs_output(out_walrus_val)
    assert [r[0] for r in rows_walrus_val] == ["r1", "r3"]

    out_walrus_target = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "pairtools",
            "select",
            "(mapq1 := 1) == 1",
            str(inp),
        ],
    ).decode("ascii")
    _, rows_walrus_target = _parse_pairs_output(out_walrus_target)
    assert [r[0] for r in rows_walrus_target] == ["r1", "r2", "r3"]
