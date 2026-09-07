# Changelog

## 3.1.0

### Added

- `exorcise-database --remove` removes experiments from a database, restoring
  the capability that `crispr-screen-viewer remove` provided before 3.0.0. The
  underlying functions were carried over in 3.0.0 but nothing reached them.

  Targets may be given either as analysis workbooks, in which case the
  experiment ID is read from the workbook's "Experiment name", or as experiment
  IDs directly. Workbooks are usually what you have to hand, since that is what
  the experiment was added with; IDs are what `experiments_metadata.csv.gz`
  lists, and what the old command took.

  ```bash
  exorcise-database --remove --out-dir db/ my-screen.xlsx
  exorcise-database --remove --out-dir db/ my-experiment-id
  ```

  Removal asks for confirmation, subject to the same rules as replacing a
  database with `--new-db`: the prompt only appears when stdin is a terminal,
  and `--force-overwrite` skips it. Without either, it refuses rather than
  proceeding. Experiments named but absent from the database are reported and
  skipped rather than being treated as an error.

  Statistics rows and both metadata tables are pruned. Gene records are left
  alone, since they are shared between experiments and hold no per-experiment
  data.

- `crispr-screen-viewer remove DB_DIR EXP_ID...` works again, with a deprecation
  warning, and is translated to
  `exorcise-database --remove --out-dir DB_DIR EXP_ID...`. The old form took the
  database directory as its first positional argument, so this needed argument
  rearrangement rather than the argument-stripping the other deprecated names
  use; `deprecation.py` now supports per-subcommand translation.

- `crispr-screen-viewer launch`, `genes` and `test` now fail with a message
  saying they went with the Dash web viewer in 3.0.0, rather than being passed
  through to `exorcise-database` as if they were arguments.

- `read_experiment_id()` reads an experiment's name from a workbook's
  "Experiment details" sheet alone. Removing an experiment needs nothing else,
  and should not fail because some sheet it will never consult is malformed.

### Changed

- The project version is now written down once, in the `VERSION` file at the
  root. Releasing means editing that file and nothing else. It had been repeated
  in four places — `bin/exorcise.R`, `bin/exorcise-trace.R`, `py/pyproject.toml`
  and `py/exorcise_screens/__init__.py` — of which only `pyproject.toml` was
  load-bearing for the image tag, so the others could drift unnoticed until
  `exorcise --version` disagreed with the tag it shipped under. They had already
  drifted once: `bin/exorcise.R` reached 3.0.1 while everything else said 3.0.0.

  Both R entry points read it through `R/version.R`, the Python package through
  `exorcise_screens.__version__`, `pyproject.toml` derives the package version
  from that attribute, and `docker/build-and-push.sh` reads the file directly.

- `bin/exorcise.R` locates its installation root rather than its `R/` directory,
  since it now needs both the modules and the version from there. `exorcise-trace`
  gained the same lookup; it sources only `R/version.R`, which depends on nothing,
  so it still does not load the Bioconductor stack.

### Fixed

- `--remove` together with `--new-db` is rejected, rather than silently letting
  one win.

## 3.0.0

crispr_tools and crispr_screen_viewer are no longer pulled from GitHub at build
time. Both have been absorbed into this repository as a single Python package,
`exorcise_screens`, holding only the functionality the toolchain actually
invokes. The Docker image now builds from the working tree.

### Commands renamed

Every old name still works and does the same thing. Each prints a warning naming
its replacement.

| Was                             | Now                 |
|---------------------------------|---------------------|
| `exorcise`                      | `exorcise`          |
| `ntByCycle`                     | `exorcise-trace`    |
| `count_reads`                   | `exorcise-count`    |
| `crispr_pipeline`               | `exorcise-analyse`  |
| `crispr-screen-viewer database` | `exorcise-database` |

The mapping is declared once, in `py/exorcise_screens/deprecation.py`, and is
used by the CLI, the shell wrappers and the Docker entrypoint alike.

### Removed

- The Dash web application: `launch`, `screen_explorer`, `multiscreen_gene_viewer`,
  `comparison_maker`, `selector_tables`, `shared_components`, `legal` and the
  `assets/` tree. With them go the `dash`, `dash-bootstrap-components`,
  `dash-bio`, `plotly` and `flask` dependencies. The database builder only ever
  needed eight small helpers out of `functions_etc.py`, none of which touch
  Dash; a single top-level `from dash import html` was dragging the entire web
  stack into every `crispr-screen-viewer database` run.
- The `launch`, `genes`, `remove` and `test` subcommands of
  `crispr-screen-viewer`. Only `database` is carried over.
- Modules unreachable from any entry point: `jacks_tools`, `plotting`, `qc`,
  `sample_down`, `pipeline_to_viewer`, `exp_class`, `tests`, `tsts`,
  `update_gene_table`, `dataset.DataSet`.
- JACKS support. It was already behind a `try/except ImportError` with a stub
  that raised, and `call_jacks` had no reachable caller.
- About 7,000 of the roughly 11,800 Python lines, none of it reachable from
  `count_reads`, `crispr_pipeline` or `crispr-screen-viewer database`.
- `fuzzyTwoLists`, a committed Go binary, and `fuzzyTwoLists.go`. See below.
- Dependencies no longer needed: `statsmodels` and `xlsxwriter` (used only by
  dead code), `scikit-learn` (declared but never imported anywhere, in either
  package), `pyaml` (only reached through `exp_class`, itself dead), and the
  `plotting` and `depreciated` extras. The latter declared `yaml`, which is not
  the PyPI name for PyYAML, so it could never have installed.

### Fixed

- `count_reads --allow-mismatch` could not work in the Docker image. It shelled
  out to `fuzzyTwoLists`, which was committed as a Mach-O x86_64 macOS build and
  so would not execute on Linux or Apple silicon. It has been reimplemented in
  Python: each library guide is expanded into the set of sequences one
  substitution away, so matching is a dictionary lookup rather than an
  all-against-all comparison. The old version also wrote three temporary files
  into the working directory and never removed them, so two concurrent runs in
  one directory corrupted each other.
- `call_mageck` built its command as a string and ran it with `shell=True`, with
  no quoting. A sample name or path containing a space or a shell metacharacter
  would break the command or inject into it. The command is now an argument list.
  The reason `shell=True` appeared necessary was that keyword arguments were
  being assembled as single strings containing a space, so each arrived as one
  argv element; they are now split properly.
- MAGeCK's exit status was ignored, so a failed comparison looked like a
  successful one until the tabulation step failed to find its output.
- `run_analyses` required `days_grown` and `cell_line_hash` positionally, but
  only `AnalysisWorkbook` produces them, so any JSON configuration raised a
  `TypeError`. They now default to empty, which is what a non-Chronos run wants.
- Result tabulation ran inside the per-analysis loop and used whichever
  `ctrl_map` had last leaked out of the inner loop over control groups. It now
  runs once, after every analysis, against the control map each analysis was
  actually given.
- `--skip-method` was matched with `in` against the raw string, so
  `--skip-method z` would silently also skip `drugz`. It is now split on commas
  like `--run-groups` and `--run-analyses` already were.
- `count_reads --debug` was parsed but never passed on, so it did nothing.
  `--quiet` claimed to imply `--just-go` but did not, so a quiet run still
  blocked on a confirmation prompt. That prompt is gone: counting no longer
  waits for confirmation, and `--just-go` is accepted and ignored.
- `count_reads -p out` (a prefix with no directory component) crashed with a
  `FileNotFoundError`, because `os.path.dirname("out")` is `""` and `os.makedirs`
  was called on it unguarded.
- Reads were identified in a FASTQ by testing each line against `[ATCGN]`, which
  a quality line can also satisfy. Records are now parsed properly, four lines
  at a time.
- `count_reads --file_type` was passed as far as the file filter but never
  reached the parser, so it had no effect on how files were read. It is now
  `--file-type`, consistent with the other flags, and is honoured.
- `tabulate_drugz` raised an unhelpful `UnboundLocalError` when no result files
  were found; it now reports which prefix it looked under, as
  `tabulate_mageck` already did.
- The database builder read `experiment_details.DOI` as an attribute, so a
  workbook without a DOI row failed with `AttributeError` rather than a message.
  A missing `Organism` behaved similarly.
- `update_database` gave a confusing error when pointed at a directory that was
  not a database. It now says which files are missing and suggests `--new-db`.
- `--verbosity 3` raised `IndexError`; it is now clamped, with a warning.
- `exorcise-database --new-db` over an existing database died with an `EOFError`
  from its own confirmation prompt whenever stdin was not a terminal, which is
  the case for `docker run` without `-it`. It now refuses to replace the database
  and says how to proceed, rather than crashing; the prompt is only shown when
  there is somebody to answer it. Ctrl+C and Ctrl+D at the prompt both cancel
  cleanly, leaving the database untouched.
- A bad `--counts-dir`, a missing counts file or a malformed workbook produced a
  traceback from `get_paths`. These are mistakes in the command, so they now
  report the problem and exit 2. Ctrl+C anywhere in any of the three commands
  exits 130 rather than printing a traceback.
- Reading a workbook emitted a pandas `FutureWarning` about silent downcasting on
  every run. Sheets are read with `dtype=object` deliberately, so the fill now
  opts out of downcasting explicitly.
- The circular import between `data_classes` and `crispr_pipeline`
  (`workbook_to_dict` imported `clargs_to_dict` from the pipeline, which imported
  the workbook) is gone: `clargs_to_dict` now lives in `util`.
- Importing anything at all from `crispr_tools` executed a package `__init__`
  that imported `count_reads`, `tools` and `exp_class`, which between them pulled
  in scipy, h5py, statsmodels, xlsxwriter, attrdictionary, PyYAML and, by way of
  `exp_class`, the pipeline itself. Nothing is imported eagerly now, so
  `exorcise-count --help` no longer waits for TensorFlow.

### Changed

- The Docker image builds from the working tree with `COPY`, not `git clone`.
  Because the build context is now the repository root, the Dockerfile must be
  named explicitly: `docker build -f docker/Dockerfile -t <tag> .`. A
  `.dockerignore` keeps the genome files and the manuscript out of the context.
- Three conda environments become two, one per language. MAGeCK and RRA now sit
  in the same environment as the code that calls them, and that environment is on
  the image's `PATH`, so the build no longer rewrites three source files with
  `sed` to hard-code absolute paths into them. The `sed` targeting MAGeCK's
  `crisprFunction.py` also hard-coded `python3.10` in a site-packages path.
- `exorcise-count` writes its per-sample counts as before but no longer prompts
  for confirmation, which made it unusable in a non-interactive container.
- Analysis and workbook validation collect every problem before reporting,
  rather than stopping at the first.
- `DETAILSTEMPLATE.xlsx` is now `example/example-analysis-workbook-template.xlsx`.
- `Hart2017_TableS2_core_genes.txt` is not carried over; nothing read it.

## 2.1.0

Refactoring release. No change to the reannotation algorithm; the pipeline
stages, their intermediate file names and the meaning of the `exo_` output
columns are all unchanged, so existing checkpoints are still valid.

### Structure

- Library code moved from `bin/` to `R/`; `bin/` now holds only the two
  executables and their wrappers.
- `bin/exorcise.R` finds its modules relative to itself, or via `EXORCISE_HOME`.
  It previously hard coded `/exorcise/bin`, so only the Docker image worked.
- BLAT and twoBitToFa are located through `EXORCISE_BLAT` and
  `EXORCISE_TWOBITTOFA`. The Dockerfile no longer rewrites the source with
  `sed` at build time.
- `psl.R` renamed to `blat.R`, `reannotate.R` to `pipeline.R` and `validate.R`
  to `options.R`.
- `data.table::fread` is no longer shadowed by a global override; reading is
  done through an explicit `read_table()`.
- Removed roughly 150 lines of commented-out code for the harmonisation and
  post-hoc features that were withdrawn in 2.0.

### Fixed

- `--exprcutoff` was compared against `class(x) == "numeric"` while always
  arriving from the command line as a character string, so any value the user
  supplied was rejected and silently replaced with the default of 10.
- Base editor `stoploss` was never reported: the consequence test compared the
  amino acid against `be_original`, a nucleotide string, instead of
  `be_original_in_aa`.
- Custom base editor modes logged `be_window_end`, which is never set, instead
  of `be_window_to`.
- `--expression` given without `--harm` aborted on a missing column part way
  through the run. It is now rejected during validation.
- In CRISPRi/a and proximity modes, padding a gene near the start of a contig
  produced a negative coordinate that was then parsed as a malformed range.
  Padding is now clamped at 1.
- Labelling non-targeting guides used `1:n`, which yields `c(1, 0)` when there
  are none and aborted on the length mismatch.
- `exo_id` was rebuilt from sequence, target and symbol after de-duplication,
  which could reintroduce identical rows. Input rows are now de-duplicated
  before keys are assigned, and the output is de-duplicated after.
- The `NA`/`""` to `"X"` sentinel round trip was applied to every column, so a
  legitimate value of `X` anywhere in the input came back as `NA` and an empty
  string became `NA`. The sentinel is gone; only `exo_symbol` is relabelled.
- The guide FASTA was written by `fwrite`, which prefixed it with a `V1` header
  line. It is now written with `writeLines`.
- Non-targeting guides were named `exo_Non-targeting1` or
  `exo_Non-targeting_1` depending on whether `--control` was given. They are
  now always `exo_Non-targeting_1`.
- Exit statuses from BLAT and twoBitToFa are checked, and a mismatch between the
  number of sequences returned and the number requested is now an error rather
  than a silent misalignment of coordinates and sequences.
- Base editor mode aborted on any codon containing an ambiguity code. Fuzzy
  codons are now resolved where the amino acid is unambiguous.
- ntByCycle: the plot title was built but never attached to the plot; a
  `--start` beyond the read length produced a descending cycle sequence; reads
  were identified by testing for `[ATCGN]`, which can also match a quality
  line; and running on more than one file blocked on a `readline` prompt that
  cannot be answered when the process is not interactive.

### Changed

- `--library` and `--priorities`, no-ops since 2.0, have been removed. Commands
  still passing them will now fail rather than silently ignore them.
- Base editor output column `be_position_in_grch38` renamed to
  `be_position_in_genome`, which is true of assemblies other than GRCh38.
- Rows whose sequence is empty or contains characters other than `[ACGTN]` are
  dropped with a warning, rather than being sent to BLAT to fail there.
- Validation reports every problem it finds at once, rather than the first.
- `--mode` accepts its argument case-insensitively and rejects a custom base
  editor window whose start is after its end.
- Bystander edit enumeration and the amino acid comparison no longer scan the
  whole table once per row.
- ntByCycle writes the counts alongside the plots as a TSV.
- Dependencies: `Biostrings` and `rtracklayer` are declared explicitly rather
  than relied upon transitively. `foreach`, `readxl`, `scales`, `stringi` and
  `plyranges` are no longer used and have been dropped. Two pins were wrong and
  would have failed to solve: `r-tidyr` was `.3.2`, which is not a version, and
  `bioconductor-rtracklayer` was `1.70.0`, which was never released.
- Two dependencies are needed without being named anywhere in the source, and are
  declared with a comment saying so, because both had been dropped as unused:
  - `R.utils`, which `data.table::fread` reaches for to decompress gzipped input.
    Most exorcise inputs are gzipped. It is no longer `library()`ed, only
    installed. `exorcise` now checks for it at startup rather than letting
    `fread` fail obscurely part way through a run.
  - `pandoc`, which `htmlwidgets::saveWidget` shells out to in order to inline
    JavaScript into a self-contained HTML file. Without it `exorcise-trace` wrote
    no interactive plot. It now falls back to a plot with its JavaScript in a
    sibling directory, and says what to install for a single file.

## 2.0.1

- Made the expression file optional; removed support for globbed inputs.

## 2.0

- Restructured; added the expression mask; removed support for harmonisation
  and post-hoc mode.

## 1.6

- Supported bystander edits in base editor mode.

## 1.5.3

- Fixed guide identifiers in base editor mode.

## 1.5.2

- Added custom base editor mode.

## 1.5.1

- Enabled inheriting arbitrary columns from the exome.

## 1.5

- Added base editor mode.

## 1.4.3

- Added arbitrary CRISPR chemistry mode with a user-supplied effect range.

## 1.4.2

- Relaxed the CRISPRi/a distance restraint; ignored genomic hits on contigs
  absent from the exome.

## 1.4

- Sped up harmonisation by vectorising it.

## 1.3

- Added CRISPRi/a support.

## 1.2

- Enforced stricter exome column naming, following the UCSC Table Browser
  format.

## 1.1

- Added `exo_id_harm` to handle non-unique input sequences.

## 1.0.2.4

- Supported exome files with comment headers.

## 1.0.2

- Warned on a potentially incorrect checkpointed PSL file.

## 1.0.1

- Improved input and output handling.

## 1.0

- First tested full release.

## 0.9

- First full release, for evaluation.
