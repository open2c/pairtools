import bioframe


def is_compatible_viewframe(
    view_df, chromsizes, check_sorting=False, raise_errors=False
):
    """
    Check if view_df is a viewframe and if
    it is compatible with the provided chromsizes.
    Parameters
    ----------
    view_df :  DataFrame
        view_df DataFrame to be validated
    chromsizes : DataFrame
        Chromsizes object to use for verification
    check_sorting : bool
        Check is regions in view_df are sorted as in
        chromosomes in cooler.
    raise_errors : bool
        raise expection instead of returning False
    Returns
    -------
    is_compatible_viewframe : bool
        True when view_df is compatible, False otherwise
    """
    try:
        try:
            _ = bioframe.is_viewframe(view_df, raise_errors=True)
        except Exception as error_not_viewframe:
            try:
                _ = bioframe.make_viewframe(view_df)
            except Exception as error_cannot_make_viewframe:
                # view_df is not viewframe and cannot be easily converted
                raise ValueError(
                    "view_df is not a valid viewframe and cannot be recovered"
                ) from error_cannot_make_viewframe
            else:
                # view_df is not viewframe, but can be converted - formatting issue ? name-column ?
                raise ValueError(
                    "view_df is not a valid viewframe, apply bioframe.make_viewframe to convert"
                ) from error_not_viewframe

        # is view_df contained inside chromsizes ?
        chromsizes = bioframe.make_viewframe(chromsizes)
        if not bioframe.is_contained(view_df, chromsizes, raise_errors=False):
            raise ValueError("View table is out of the bounds of chromosizes.")

        # is view_df sorted by coord and chrom order as in cooler ?
        if check_sorting:
            if not bioframe.is_sorted(view_df, chromsizes, df_view_col="chrom"):
                raise ValueError(
                    "regions in the view_df must be sorted by coordinate"
                    " and chromosomes order as as in the chromsizes."
                )

    except Exception as e:
        if raise_errors:
            raise ValueError("view_df is not compatible, or not a viewframe") from e
        else:
            # something went wrong: it's not a viewframe
            return False
    else:
        # no exceptions were raised: it's a compatible viewframe
        return True


def read_viewframe_from_file(
    view_fname,
    verify_chromsizes=None,
    check_sorting=False,
):
    """
    Read a BED file with regions that conforms
    a definition of a viewframe (non-overlaping, unique names, etc).
    Parameters
    ----------
    view_fname : str
        Path to a BED file with regions.
    verify_chromsizes : DataFrame | None
        Chromsizes for bound checking
        No checks are done when None.
    check_sorting : bool
        Check is regions in view_df are sorted as in
        chromosomes in cooler.
    Returns
    -------
    view_df : pd.DataFrame
        DataFrame with the viewframe
    """

    # read BED file assuming bed4/3 formats (with names-columns and without):
    try:
        view_df = bioframe.read_table(view_fname, schema="bed4", index_col=False)
    except Exception as err_bed4:
        try:
            view_df = bioframe.read_table(view_fname, schema="bed3", index_col=False)
        except Exception as err_bed3:
            raise ValueError(
                f"{view_fname} is not a BED file with 3 or 4 columns"
            ) from err_bed4

    # Convert view dataframe to viewframe:
    try:
        view_df = bioframe.make_viewframe(view_df)
    except ValueError as e:
        raise ValueError(
            "View table is incorrect, please, comply with the format. "
        ) from e

    if verify_chromsizes is not None:
        try:
            _ = is_compatible_viewframe(
                view_df, verify_chromsizes, check_sorting, raise_errors=True
            )
        except Exception as e:
            raise ValueError("view_df is not compatible with the cooler") from e
        else:
            # view_df is compaible, returning
            return view_df
    else:
        # no cooler for checking, returning
        return view_df
