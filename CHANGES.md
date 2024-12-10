### 1.1.1 (2024-12-10) ###

Bugfixes:
- Migrating to pyproject.toml + cibuildwheel. pairtools will now release binary wheels for Linux. --no-build-isolation is a mandatory flag now.
- Require Cython during build to avoid the "circular import" bug.
- fix API incomplete functionality for read-side detection by @agalitsyna 

**Full Changelog**: https://github.com/open2c/pairtools/compare/v1.1.0...v1.1.1

### 1.1.0 (2024-04-23) ###
Major bugfixes:
- Fix a major bug in sort that previously broke the sorting order. This bug was introduced in recent versions of pairtools #230
- Fix a major bug in dedup that caused pair duplication and broken sorting order in non-Cython backends

New features:
- stats: calculate the distance of P(s) divergence between pairs of different directionalities #222
- dedup: allow column names in all backends, and allow sorting by arbitrary columns #162

New behavior and default settings:
- dedup: turn mark-dups on by default #211
- parse: change the default --walks-policy to 5unique
- parse: pair types are now always in upper case. Previously, letters in pair types were converted to lowercase if the corresponding side contained chimeric alignments.

Minor bugfixes:
- dedup: allow inputs with quotes #194
- dedup: allow empty input pairs file #201
- stats: minor bugfixes #200

Documentation:
- a new notebook with the statistics of distances between PCR duplicates #233
- clean up phase walkthrough #218
- a new chapter on building workflows with pairtools #219 #226 #231
- a major cleanup

Code updates:
- make pairsio.py to read (and, in the future, write) .pairs files #195
- minor refactoring of parse #223

New Contributors:
- @hkariti made their first contribution in #194

### 1.0.3 (2023-11-20) ###
- [x] `pairtools dedup`: update default chunksize to 10,000 to prevent memory overflow on datasets with high duplication rate 

### 1.0.2 (2022-11-XX) ###

- [x] `pairtools select` regex update 
(string substitutions failed when the column name was a substring of another)

- [x] Warnings capture in dedup: pairs lines are always split after rstrip newline

- [x] Important fixes of splitting schema

- [x] Dedup comment removed (failed when the read qualities contained "#")

- [x] Remove dbist build out of wheel

- [x] pairtools scaling: fixed an issue with scaling maximum range value https://github.com/open2c/pairtools/issues/150#issue-1439106031 

### 1.0.1 (2022-09-XX) ###

- [x] Fixed issue with pysam dependencies on pip and conda

- [x] pytest test engine instead of nose

- [x] Small fixes in teh docs and scaling

### 1.0.0 (2022-08-XX) ###

This is a major release of pairtools since last release (April 2019!)

#### Post merge:

- [x] sphinx docs update with incorporated walkthroughs

#### New tools:
- [x] parse2 module with CLI for parsing complex walks
- [x] scaling and header modules with CLI

#### Fixes by modules:

pairtools dedup
- [x] finalize detection of optical duplicates https://github.com/open2c/pairtools/issues/106 and https://github.com/open2c/pairtools/issues/59, also related to  https://github.com/open2c/pairtools/issues/54 
- [x] chunked dedup by @Phlya 
- [x] improvement of dedup to include reporting of the parent readID by @Phlya and @agalitsyna

pairtools stats/scaling
- [x] split dedup stats and regular stats
- [x] output chromosome size to the stats output https://github.com/open2c/pairtools/issues/83 
- [x] pairtools stats: YAML output? https://github.com/open2c/pairtools/issues/111  and https://github.com/open2c/pairtools/issues/79
- [x] pairtools scaling tool which takes into account chromosome sizes: https://github.com/open2c/pairtools/issues/81,  https://github.com/open2c/pairtools/issues/56? 

pairtools parse
- [x] parse complex walks engine and tools: https://github.com/open2c/pairtools/pull/109
- [x] stdin and stdout reporting defaults: https://github.com/open2c/pairtools/issues/48 
- [x] flipping issue: https://github.com/open2c/pairtools/issues/91 

pairtools phase
- [x] make work with both pip and github versions of bwa: https://github.com/open2c/pairtools/pull/114

pairtools restrict
- [x] Handle empty pairs with "!" chromosomes: https://github.com/open2c/pairtools/issues/76 
- [x] Problem with restriction sites header/first rfrag: https://github.com/open2c/pairtools/issues/73 
- [x] Suggestions by @golobor: https://github.com/open2c/pairtools/issues/16

pairtools merge
- [x] do not require sorting? https://github.com/open2c/pairtools/issues/23 
- [x] headers handling: https://github.com/open2c/pairtools/issues/18

#### General improvements:

Headers maintenance
- [x] allow adding a header to a headerless file https://github.com/open2c/pairtools/issues/119
or broader addition of the headed module, draft: https://github.com/open2c/pairtools/pull/121 

Code maintenance
- [x] transfer pairlib into sandbox of pairtools lib
- [x] separate cli and lib
- [x] Remove OrderedDict: https://github.com/open2c/pairtools/issues/113 
- [x] Clean up deprecation warnings, e.g. https://github.com/open2c/pairtools/issues/71
- [x] Fix input errors without explanations, e.g. https://github.com/open2c/pairtools/issues/61 

#### Specific changes: 

Docs improvements
- [x] pairtools walkthrough
- [x] phasing walkthrough
- [x] parse docs update

Tests proposals
- [x] add tests for dedup @Phlya : https://github.com/open2c/pairtools/issues/5
- [x] add tests for stats, and merge: https://github.com/open2c/pairtools/issues/5

Enhancements
- [x] add summaries: https://github.com/open2c/pairtools/pull/105 
- [x] support of [bwa mem2]( https://github.com/bwa-mem2/bwa-mem2), which is 2-3 times faster than usual bwa mem: https://github.com/open2c/pairtools/discussions/118
- [x] I/O single utility instead of repetitive code in each module


### 0.3.1 (2021-02-XX) ###

* sample: a new tool to select a random subset of pairs
* parse: add --readid-transform to edit readID
* parse: add experimental --walk-policy all (note: it will be moved 
  to a separate tool in future!) 
* all tools: use bgzip if pbgzip not available

Internal changes:
* parse: move most code to a separate _parse module
* _headerops: add extract_chromosomes(header)  
* all tools: drop py3.5 support
* switch from travis CI to github actions

### 0.3.0 (2019-04-23) ###

* parse: tag pairs with missing FASTQ/SAM on one side as corrupt, pair type "XX"

### 0.2.2 (2019-01-07) ###

* sort: enable lz4c compression of sorted chunks by default

### 0.2.1 (2018-12-21) ###

* automatically convert mapq1 and mapq2 to int in `select`

### 0.2.0 (2018-09-03) ###

* add the `flip` tool

### 0.1.1 (2018-07-19) ###

* Bugfix: include _dedup.pyx in the Python package

### 0.1.0 (2018-07-19) ###

* First release.
