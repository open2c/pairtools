from ..lib import fileio, pairsam_format, headerops
import re, fnmatch

# Create environment of important functions:
wildcard_library = {}


def wildcard_match(x, wildcard):
    if wildcard not in wildcard_library:
        regex = fnmatch.translate(wildcard)
        reobj = re.compile(regex)
        wildcard_library[wildcard] = reobj
    return wildcard_library[wildcard].fullmatch(x)


csv_library = {}


def csv_match(x, csv):
    if csv not in csv_library:
        csv_library[csv] = set(csv.split(","))
    return x in csv_library[csv]


regex_library = {}


def regex_match(x, regex):
    if regex not in regex_library:
        reobj = re.compile(regex)
        regex_library[regex] = reobj
    return regex_library[regex].fullmatch(x)


# Define default data types:
TYPES = {"pos1": "int", "pos2": "int", "mapq1": "int", "mapq2": "int"}


def evaluate_stream(
    headerless_stream, condition, column_names, type_cast=(), startup_code=None
):
    """
    Evaluate expression for the input headerless stream.

    Parameters
    ----------
    headerless_stream
    condition
    type_cast: Cast a given column to a given type. By default, only pos and mapq
                are cast to int, other columns are kept as str. Type: tupe of two strings.
    startup_code: An auxiliary code to execute before filtering.
                Use to define functions that can be evaluated in the CONDITION statement

    ========
    Writes the output to one of two streams (regular or rest)

    """

    # Define data types:
    TYPES.update(dict(type_cast))

    # Execute startup code:
    if startup_code is not None:
        exec(startup_code, globals())

    for i, col in enumerate(column_names):
        if col in TYPES:
            col_type = TYPES[col]
            condition = re.sub(r"\b%s\b" % col , "{}(COLS[{}])".format(col_type, i), condition)
            #condition.replace(col, "{}(COLS[{}])".format(col_type, i))
        else:
            condition = re.sub(r"\b%s\b" % col, "COLS[{}]".format(i), condition)
            #condition = condition.replace(col, "COLS[{}]".format(i))

    # Compile the filtering expression:
    match_func = compile(condition, "<string>", "eval")

    for line in headerless_stream:
        COLS = line.rstrip().split(pairsam_format.PAIRSAM_SEP)

        # Evaluate filtering expression:
        filter_passed = eval(match_func)

        # Produce the output:
        yield filter_passed, line


def evaluate_df(df, condition, type_cast=(), startup_code=None, engine="pandas"):
    """
    Evaluate expression for the input headerless stream.

    Parameters
    ----------
    df: input dataframe for evaluation
    condition: condition to evaluate
    type_cast: additional types transformations, if different from default
    startup_code: An auxiliary code to execute before filtering.
                Use to define functions that can be evaluated in the CONDITION statement

    ========
    Writes the output to one of two streams (regular or rest)

    """

    # Define data types:
    TYPES.update(dict(type_cast))

    # Execute startup code:
    if startup_code is not None:
        exec(startup_code, globals())

    # Set up the column formats:
    for col in df.columns:
        if col in TYPES.keys():
            if not str(df.dtypes[col]) != TYPES[col]:
                df[col] = df[col].astype(TYPES[col])

    if engine == "pandas":
        try:
            filter_passed_output = df.eval(condition)
        except ValueError as e:
            raise ValueError(f"Try passing engine python to fix this: {e}")
    else:
        # Set up the columns indexing
        for i, col in enumerate(df.columns):
            condition = re.sub(r"\b%s\b" % col, "COLS[{}]".format(i), condition)
            #condition = condition.replace(col, "COLS[{}]".format(i))

        filter_passed_output = []
        match_func = compile(condition, "<string>", "eval")
        for i, r in df.iterrows():
            COLS = r.values

            # Evaluate filtering expression:
            filter_passed = eval(match_func)
            filter_passed_output.append(True if filter_passed else False)

    return filter_passed_output
