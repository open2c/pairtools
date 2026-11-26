from ..lib import fileio, pairsam_format, headerops
import ast
import copy
import keyword
import re, fnmatch
import numpy as np

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


def region_match(chrom, pos, region_chrom, region_start=-1, region_end=-1):
    if region_end == -1:
        region_end = np.inf
    return chrom == region_chrom and region_start <= pos <= region_end


# Define default data types:
TYPES = {"pos1": "int", "pos2": "int", "mapq1": "int", "mapq2": "int"}


def _extract_target_names(target):
    """Return variable names bound (Store ctx) by a comprehension target."""
    if isinstance(target, ast.Name):
        return {target.id} if isinstance(target.ctx, ast.Store) else set()
    if isinstance(target, (ast.Tuple, ast.List)):
        names = set()
        for elt in target.elts:
            names.update(_extract_target_names(elt))
        return names
    if isinstance(target, ast.Starred):
        return _extract_target_names(target.value)
    return set()


class _ColumnNameBinder(ast.NodeTransformer):
    """AST transformer that binds column names while respecting Python scoping."""

    def __init__(self, replacements):
        self._replacements = replacements
        self._parsed_cache = {}
        self._shadowed = set()

    def _get_replacement_node(self, col):
        node = self._parsed_cache.get(col)
        if node is None:
            node = ast.parse(self._replacements[col], mode="eval").body
            self._parsed_cache[col] = node
        return copy.deepcopy(node)

    def visit_Name(self, node):
        if (
            isinstance(node.ctx, ast.Load)
            and node.id in self._replacements
            and node.id not in self._shadowed
        ):
            return ast.copy_location(self._get_replacement_node(node.id), node)
        return node

    def visit_Lambda(self, node):
        # Default argument expressions are evaluated in the enclosing scope.
        node.args.defaults = [self.visit(d) for d in node.args.defaults]
        node.args.kw_defaults = [
            self.visit(d) if d is not None else None for d in node.args.kw_defaults
        ]

        param_names = {
            arg.arg
            for arg in (
                list(node.args.posonlyargs)
                + list(node.args.args)
                + list(node.args.kwonlyargs)
            )
        }
        if node.args.vararg is not None:
            param_names.add(node.args.vararg.arg)
        if node.args.kwarg is not None:
            param_names.add(node.args.kwarg.arg)

        prev_shadowed = self._shadowed
        self._shadowed = prev_shadowed | param_names
        try:
            node.body = self.visit(node.body)
        finally:
            self._shadowed = prev_shadowed
        return node

    def _visit_comprehension(self, node, body_attrs):
        # The outermost generator's iter is evaluated in the enclosing scope;
        # all comprehension targets are local to the comprehension scope.
        if node.generators:
            node.generators[0].iter = self.visit(node.generators[0].iter)

        target_names = set()
        for gen in node.generators:
            target_names.update(_extract_target_names(gen.target))

        prev_shadowed = self._shadowed
        self._shadowed = prev_shadowed | target_names
        try:
            for idx, gen in enumerate(node.generators):
                if idx > 0:
                    gen.iter = self.visit(gen.iter)
                gen.target = self.visit(gen.target)
                gen.ifs = [self.visit(if_clause) for if_clause in gen.ifs]
            for attr in body_attrs:
                setattr(node, attr, self.visit(getattr(node, attr)))
        finally:
            self._shadowed = prev_shadowed
        return node

    def visit_ListComp(self, node):
        return self._visit_comprehension(node, ("elt",))

    def visit_SetComp(self, node):
        return self._visit_comprehension(node, ("elt",))

    def visit_GeneratorExp(self, node):
        return self._visit_comprehension(node, ("elt",))

    def visit_DictComp(self, node):
        return self._visit_comprehension(node, ("key", "value"))


def _compile_condition(condition, replacements):
    """Parse, bind column names in the AST, and compile a Python condition.

    Replaces only ``ast.Name`` nodes in ``ast.Load`` context that are not
    shadowed by enclosing lambda parameters or comprehension targets.
    Quoted literals, f-string debug prefixes (``f"{col=}"``), attribute names
    (``obj.col``), keyword argument names (``fn(col=...)``), and walrus
    targets (``(col := ...)``) are preserved by Python's AST semantics across
    Python 3.9+.
    """
    tree = ast.parse(condition, mode="eval")
    if replacements:
        tree = _ColumnNameBinder(replacements).visit(tree)
        ast.fix_missing_locations(tree)
    return compile(tree, "<string>", "eval")


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

    # Define data types without mutating module-level defaults:
    _types = TYPES.copy()
    _types.update(dict(type_cast))

    # Execute startup code:
    if startup_code is not None:
        exec(startup_code, globals())

    replacements = {}
    for i, col in enumerate(column_names):
        if (
            col not in replacements
            and isinstance(col, str)
            and col.isidentifier()
            and not keyword.iskeyword(col)
        ):
            if col in _types:
                replacements[col] = "{}(COLS[{}])".format(_types[col], i)
            else:
                replacements[col] = "COLS[{}]".format(i)
    # Compile the filtering expression:
    match_func = _compile_condition(condition, replacements)

    for line in headerless_stream:
        COLS = line.rstrip("\n").split(pairsam_format.PAIRSAM_SEP)

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

    # Define data types without mutating module-level defaults:
    _types = TYPES.copy()
    _types.update(dict(type_cast))

    # Execute startup code:
    if startup_code is not None:
        exec(startup_code, globals())

    # Set up the column formats:
    for col in df.columns:
        if col in _types.keys():
            if not str(df.dtypes[col]) != _types[col]:
                df[col] = df[col].astype(_types[col])

    if engine == "pandas":
        try:
            filter_passed_output = df.eval(condition)
        except ValueError as e:
            raise ValueError(f"Try passing engine python to fix this: {e}")
    else:
        # Set up the columns indexing
        replacements = {}
        for i, col in enumerate(df.columns):
            if (
                col not in replacements
                and isinstance(col, str)
                and col.isidentifier()
                and not keyword.iskeyword(col)
            ):
                replacements[col] = "COLS[{}]".format(i)

        filter_passed_output = []
        match_func = _compile_condition(condition, replacements)
        for i, r in df.iterrows():
            COLS = r.values

            # Evaluate filtering expression:
            filter_passed = eval(match_func)
            filter_passed_output.append(True if filter_passed else False)

    return filter_passed_output
