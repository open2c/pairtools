import pandas as pd

from . import fileio, headerops

def read_pairs(pairs, nproc=3, cmd_in=None, **kwargs):
    """
    Reads a file with .pairs format and returns a header, a dataframe of pairs, and chromsizes.

    Parameters:
        pairs (str or file-like object): A path to a .pairs file to read or an open file-like object/handle.
        nproc (int): Number of processes to use for reading the file. Default is 3.
        cmd_in (str): The command to be used for reading the file. Default is None.

        **kwargs: Additional keyword arguments to be passed to pd.read_csv. Useful options include:
            - chunksize (int): If specified, return an iterable object of type TextFileReader that reads in chunks of lines.
            - usecols (list-like or callable): Return a subset of the columns. If list-like, all elements must either be positional or strings. If callable, the callable function will be evaluated against the column names, returning names where the callable function evaluates to True.

    Returns:
        tuple: A tuple containing the following elements:
            - pairs_df (pd.DataFrame): A pandas DataFrame with pairs.
            - header (list of str): The original header of the pairs file.
            - chromsizes (dict): A dictionary containing chromosome sizes extracted from the header.
    """
    pairs_stream = (
        fileio.auto_open(
            pairs,
            mode="r",
            nproc=nproc,
            command=cmd_in,
        )
        if isinstance(pairs, str)
        else pairs
    )

    header, pairs_body = headerops.get_header(pairs_stream)
    cols = headerops.extract_column_names(header)

    chromsizes = headerops.extract_chromsizes(header)

    pairs_df = pd.read_csv(
        pairs_body,
        header=None,
        names=cols,
        sep="\t",
        dtype={"chrom1": str, "chrom2": str},
        **kwargs
    )

    return pairs_df, header, chromsizes