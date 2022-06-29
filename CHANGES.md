### 0.4.0 (2022-06-XX) ###

This is a major release

* parse
- [x] parse complex walks engine and tools: https://github.com/open2c/pairtools/pull/109
- [x] stdin and stdout reporting defaults: https://github.com/open2c/pairtools/issues/48 
- [x] flipping issue: https://github.com/open2c/pairtools/issues/91 

Reporting:
- pair types []
- new fields
- report all

- changes in defaults

* parse2

* pairtools dedup
- [x] finalize detection of optical duplicates https://github.com/open2c/pairtools/issues/106 and https://github.com/open2c/pairtools/issues/59, also related to  https://github.com/open2c/pairtools/issues/54 
- [x] chunked dedup by @Phlya 
- [x] improvement of dedup to include reporting of the parent readID by @Phlya and @agalitsyna

* dedup: implemented in chunks with two new backends ("scipy", "sklearn"). Now allows to
  record the readID of the retained "parent" read from a duplicate cluster in an extra
  field in the file with duplicates. New backends rely on the header to define column
  oder in the file, specification through CLI arguments works for the "cython" backend,
  but it will be removed in a future version.
  Note that with non-zero max-mismatch the behaviour of the new backends can be
  different from the old "cython": now duplication is transitive (i.e. if read A is a
  duplicate of read B, and read B - of read C, reads A and C are now considered
  duplicates).

Note: this is the header post connecting multiple issues, feel free to update and improve!


pairtools stats/scaling
- [x] split dedup stats and regular stats
- [x] output chromosome size to the stats output https://github.com/open2c/pairtools/issues/83 
- [x] pairtools stats: YAML output? https://github.com/open2c/pairtools/issues/111  and https://github.com/open2c/pairtools/issues/79
- [x] pairtools scaling tool which takes into account chromosome sizes: https://github.com/open2c/pairtools/issues/81,  https://github.com/open2c/pairtools/issues/56? 


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

#### Specific proposals: 

Docs improvements
- [x] pairtools walkthrough
- [x] phasing walkthrough
- [x] parse docs update
- [x] sphinx docs update with incorporated walkthroughs


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
