# -*- coding: utf-8 -*-
from pairtools.lib import headerops

import pytest


import pytest
from pairtools.lib.headerops import (
    canonicalize_columns,
    get_column_index,
    extract_column_names,
)

def test_canonicalize_columns():
    # Test basic canonicalization
    assert canonicalize_columns(['chr1', 'chr2']) == ['chrom1', 'chrom2']
    assert canonicalize_columns(['chrom1', 'chrom2']) == ['chrom1', 'chrom2']
    assert canonicalize_columns(['pt', 'other']) == ['pair_type', 'other']
    
    # Test mixed case
    assert canonicalize_columns(['Chr1', 'CHR2']) == ['chrom1', 'chrom2']
    assert canonicalize_columns(['CHR1', 'chr2']) == ['chrom1', 'chrom2']
    
    # Test no changes needed
    assert canonicalize_columns(['readID', 'pos1']) == ['readID', 'pos1']
    
    # Test empty input
    assert canonicalize_columns([]) == []
    
    # Test all known aliases
    assert canonicalize_columns(['chr1', 'Chr2', 'PT']) == ['chrom1', 'chrom2', 'pair_type']

def test_get_column_index():
    # Setup test columns
    columns = ['readID', 'chr1', 'pos1', 'chr2', 'pos2', 'strand1', 'strand2', 'pair_type']
    
    # Test string lookup - direct matches
    assert get_column_index(columns, 'chr1') == 1
    assert get_column_index(columns, 'pos2') == 4
    assert get_column_index(columns, 'pair_type') == 7
    
    # Test string lookup - canonicalized matches
    assert get_column_index(columns, 'chrom1') == 1
    assert get_column_index(columns, 'CHROM2') == 3
    assert get_column_index(columns, 'PT') == 7
    
    # Test case insensitive matches
    assert get_column_index(columns, 'CHR1') == 1
    assert get_column_index(columns, 'ChR2') == 3
    
    # Test integer lookup
    assert get_column_index(columns, 0) == 0
    assert get_column_index(columns, 3) == 3
    assert get_column_index(columns, 7) == 7
    
    # Test error cases
    with pytest.raises(ValueError, match="Column 'nonexistent' not found"):
        get_column_index(columns, 'nonexistent')
        
    with pytest.raises(ValueError, match="Column index 100 out of range"):
        get_column_index(columns, 100)
        
    with pytest.raises(AttributeError, match="Column spec must be string or integer"):
        get_column_index(columns, 3.14)

def test_integration_with_extract_column_names():
    # Test with actual header format
    header = [
        "## pairs format v1.0",
        "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 pair_type",
        "#chromsize: chr1 1000",
        "#chromsize: chr2 800"
    ]
    
    columns = extract_column_names(header)
    assert columns == ['readID', 'chr1', 'pos1', 'chr2', 'pos2', 'strand1', 'strand2', 'pair_type']
    
    # Test canonicalized column access
    assert get_column_index(columns, 'chrom1') == 1
    assert get_column_index(columns, 'chrom2') == 3
    assert get_column_index(columns, 'pt') == 7
    
    # Test with alternative header format
    header2 = [
        "## pairs format v1.0",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type",
    ]
    columns2 = extract_column_names(header2)
    assert get_column_index(columns2, 'chr1') == 1
    assert get_column_index(columns2, 'chr2') == 3

def test_edge_cases():
    # Test empty columns
    with pytest.raises(ValueError):
        get_column_index([], 'chrom1')
    
    # Test invalid column spec type
    with pytest.raises(AttributeError):
        get_column_index(['a', 'b', 'c'], 3.14)  # float not supported
        
    # Test negative indices
    assert get_column_index(['a', 'b', 'c'], -1) == 2  # Python-style negative indexing


def test_make_standard_header():
    header = headerops.make_standard_pairsheader()

    assert any([l.startswith("## pairs format") for l in header])
    assert any([l.startswith("#shape") for l in header])
    assert any([l.startswith("#columns") for l in header])

    header = headerops.make_standard_pairsheader(
        chromsizes=[("b", 100), ("c", 100), ("a", 100)]
    )

    assert sum([l.startswith("#chromsize") for l in header]) == 3


def test_samheaderops():
    header = headerops.make_standard_pairsheader()
    samheader = [
        "@SQ\tSN:chr1\tLN:100",
        "@SQ\tSN:chr2\tLN:100",
        "@SQ\tSN:chr3\tLN:100",
        "@PG\tID:bwa\tPN:bwa\tCL:bwa",
        "@PG\tID:bwa-2\tPN:bwa\tCL:bwa\tPP:bwa",
    ]
    header_with_sam = headerops.insert_samheader(header, samheader)

    assert len(header_with_sam) == len(header) + len(samheader)
    for l in samheader:
        assert any([l2.startswith("#samheader") and l in l2 for l2 in header_with_sam])

    # test adding new programs to the PG chain
    header_extra_pg = headerops.append_new_pg(header_with_sam, ID="test", PN="test")

    # test if all lines got transferred
    assert all([(old_l in header_extra_pg) for old_l in header_with_sam])
    # test if one PG got added
    assert len(header_extra_pg) == len(header_with_sam) + 1

    # test if the new PG has PP matching the ID of one of already existing PGs
    new_l = [l for l in header_extra_pg if l not in header_with_sam][0]
    pp = [f[3:] for f in new_l.split("\t") if f.startswith("PP:")][0]
    assert (
        len(
            [
                l
                for l in header_extra_pg
                if l.startswith("#samheader") and ("\tID:{}\t".format(pp) in l)
            ]
        )
        == 1
    )


def test_merge_pairheaders():
    headers = [["## pairs format v1.0"], ["## pairs format v1.0"]]
    merged_header = headerops._merge_pairheaders(headers)
    assert merged_header == headers[0]

    headers = [["## pairs format v1.0", "#a"], ["## pairs format v1.0", "#b"]]
    merged_header = headerops._merge_pairheaders(headers)
    assert merged_header == ["## pairs format v1.0", "#a", "#b"]

    headers = [
        ["## pairs format v1.0", "#chromsize: chr1 100", "#chromsize: chr2 200"],
        ["## pairs format v1.0", "#chromsize: chr1 100", "#chromsize: chr2 200"],
    ]
    merged_header = headerops._merge_pairheaders(headers)
    assert merged_header == headers[0]


def test_merge_different_pairheaders():
    with pytest.raises(Exception):
        headers = [["## pairs format v1.0"], ["## pairs format v1.1"]]
        merged_header = headerops._merge_pairheaders(headers)


def test_force_merge_pairheaders():
    headers = [
        ["## pairs format v1.0", "#chromsize: chr1 100"],
        ["## pairs format v1.0", "#chromsize: chr2 200"],
    ]
    merged_header = headerops._merge_pairheaders(headers, force=True)
    assert merged_header == [
        "## pairs format v1.0",
        "#chromsize: chr1 100",
        "#chromsize: chr2 200",
    ]


def test_merge_samheaders():
    headers = [
        ["@HD\tVN:1"],
        ["@HD\tVN:1"],
    ]
    merged_header = headerops._merge_samheaders(headers)
    assert merged_header == headers[0]

    headers = [
        [
            "@HD\tVN:1",
            "@SQ\tSN:chr1\tLN:100",
            "@SQ\tSN:chr2\tLN:100",
        ],
        [
            "@HD\tVN:1",
            "@SQ\tSN:chr1\tLN:100",
            "@SQ\tSN:chr2\tLN:100",
        ],
    ]
    merged_header = headerops._merge_samheaders(headers)
    assert merged_header == headers[0]

    headers = [
        [
            "@HD\tVN:1",
            "@PG\tID:bwa\tPN:bwa\tPP:cat",
        ],
        [
            "@HD\tVN:1",
            "@PG\tID:bwa\tPN:bwa\tPP:cat",
        ],
    ]
    merged_header = headerops._merge_samheaders(headers)
    print(merged_header)
    assert merged_header == [
        "@HD\tVN:1",
        "@PG\tID:bwa-1\tPN:bwa\tPP:cat-1",
        "@PG\tID:bwa-2\tPN:bwa\tPP:cat-2",
    ]


def test_merge_headers():
    headers = [
        [
            "## pairs format v1.0",
            "#samheader: @HD\tVN:1",
            "#samheader: @SQ\tSN:chr1\tLN:100",
            "#samheader: @SQ\tSN:chr2\tLN:100",
        ]
    ] * 2

    merged_header = headerops.merge_headers(headers)
    assert merged_header == headers[0]
