# -*- coding: utf-8 -*-
from pairtools.lib import headerops

import pytest


import pytest
from pairtools.lib.headerops import (
    standardize_column,
    get_column_index,
    extract_column_names,
)

def test_standardize_column():
    # Test basic standardization
    assert standardize_column('chr1') == 'chrom1'
    assert standardize_column('chr2') == 'chrom2'
    assert standardize_column('pt') == 'pair_type'
    
    # Test no changes needed
    assert standardize_column('chrom1') == 'chrom1'
    assert standardize_column('readID') == 'readID'
    
    # Test unknown columns remain unchanged
    assert standardize_column('unknown') == 'unknown'


def test_get_column_index():
    # Setup test columns
    columns = ['readID', 'chr1', 'pos1', 'chr2', 'pos2', 'strand1', 'strand2', 'pair_type']
    
    # Test string lookup - direct matches
    assert get_column_index(columns, 'chr1') == 1
    assert get_column_index(columns, 'pos2') == 4
    assert get_column_index(columns, 'pair_type') == 7
    
    # Test string lookup - standardized matches
    # Note: 'chrom1' standardizes to 'chrom1' (no change if not in alias map depending on map content, 
    # but here 'chr1' is in columns, so we search for that. 
    # Wait, 'chrom1' is NOT in columns list above. 'chr1' is.
    # standardized('chrom1') -> 'chrom1'. 'chrom1' is not in columns.
    # So searching for 'chrom1' should FAIL unless we also standardize the columns list before search,
    # which get_column_index DOES NOT do (it expects caller to do it).
    # actually, looking at the code I wrote for get_column_index:
    # it tries direct match, then "standardized" match.
    # if I pass 'chr1', standardized is 'chrom1'.
    # if I pass 'chrom1', standardized is 'chrom1'.
    # The columns list has 'chr1'. 
    # direct match 'chr1' -> found.
    # direct match 'chrom1' -> not found. standardized 'chrom1' -> 'chrom1' -> not found.
    # So 'chrom1' should fail lookup in 'columns' list ['chr1', ...].
    
    # Let's adjust the test case to match HOW it is used. 
    # In the code, we standardize the header columns FIRST.
    # So the columns list passed to get_column_index usually has STANDARDIZED names.
    
    std_columns = ['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type']
    
    assert get_column_index(std_columns, 'chrom1') == 1
    assert get_column_index(std_columns, 'chr1') == 1  # chr1 -> chrom1 (standardized) -> found
    
    assert get_column_index(std_columns, 'pt') == 7    # pt -> pair_type -> found
    assert get_column_index(std_columns, 'pair_type') == 7

    
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
    # The extraction just gets names.
    assert columns == ['readID', 'chr1', 'pos1', 'chr2', 'pos2', 'strand1', 'strand2', 'pair_type']
    
    # Standardize them
    std_columns = [standardize_column(c) for c in columns]
    assert std_columns == ['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type']
    
    # Test lookup
    assert get_column_index(std_columns, 'chrom1') == 1
    assert get_column_index(std_columns, 'chr1') == 1
    assert get_column_index(std_columns, 'pt') == 7
    
    # Test with alternative header format (already standardized names)
    header2 = [
        "## pairs format v1.0",
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type",
    ]
    columns2 = extract_column_names(header2)
    std_columns2 = [standardize_column(c) for c in columns2]
    
    assert get_column_index(std_columns2, 'chr1') == 1
    assert get_column_index(std_columns2, 'chrom1') == 1


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
