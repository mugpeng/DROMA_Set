# DROMA SQLite Benchmark Summary

Date: 2026-04-16

## Background

Observed issue:

- SQLite-backed DROMA loading is noticeably slower on Linux than on macOS.
- The main concern is large omics reads, especially `loadMolecularProfiles()` on wide tables.

Constraint for this round:

- Do not modify `DROMA_Set`.
- Build benchmark and optimization logic externally under `DROMA_benchmark`.
- Validate whether the slowdown is caused by SQLite interaction patterns rather than OS alone.

## Benchmark Setup

Real database used:

- `/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite`

Project used for real benchmark:

- `CCLE`

Key real tables:

- `CCLE_mRNA`: `18881 x 1407`
- `CCLE_meth`: `20192 x 844`
- `CCLE_proteinrppa`: `214 x 900`
- `CCLE_mutation_gene`: `779882` rows
- `CCLE_drug`: `24 x 504`

Benchmark implementation:

- External wrapper and helper code under `DROMA_benchmark/`
- No changes made to `DROMA_Set`
- Regression checks passed on the real database before timing comparison

Relevant benchmark artifacts:

- `DROMA_benchmark/droma_sqlite_wrappers.R`
- `DROMA_benchmark/droma_sqlite_helpers.R`
- `DROMA_benchmark/output/sqlite_benchmark_results.csv`

## What Was Optimized Externally

The wrapper layer focused on low-risk SQLite access improvements:

- connection reuse
- caching `dbListTables()`
- caching `PRAGMA table_info(...)`
- caching repeated `SELECT DISTINCT ...`
- selecting only needed sample columns from wide matrix tables
- avoiding `SELECT *` followed by heavy R-side column trimming when only sample subsets are needed

No public API behavior was intentionally changed in the wrapper comparison.

## Real Benchmark Results

Using `mean_sec` from the real `CCLE` benchmark:

### 1. Object creation / metadata load

- `createDromaSetFromDatabase`
  - original: `0.110 s`
  - wrapper: `0.116 s`

Interpretation:

- initial object creation is not the main bottleneck
- metadata loading overhead exists, but it is not the dominant cost in this workload

### 2. Full mRNA matrix load

- `loadMolecularProfiles(..., feature_type = "mRNA", return_data = TRUE)`
  - original: `3.048 s`
  - wrapper: `0.947 s`
  - wrapper warm: `0.775 s`

Interpretation:

- full-width omics loading is a major bottleneck
- the wrapper gives about `3.2x` improvement on the real `CCLE_mRNA` table

### 3. mRNA feature subset load

- `loadMolecularProfiles(..., feature_type = "mRNA", select_features = subset)`
  - original: `0.018 s`
  - wrapper: `0.005 s`

Interpretation:

- feature filtering already benefits from the existing `feature_id` index
- improvement exists, but this is not the largest pain point

### 4. mRNA sample subset load

- `loadMolecularProfiles(..., feature_type = "mRNA", samples = subset)`
  - original: `3.254 s`
  - wrapper: `0.054 s`
  - wrapper warm: `0.054 s`

Interpretation:

- this is the clearest bottleneck
- current logic behaves like a near-full table read even when only a sample subset is needed
- external optimized path improves this by about `60x`

### 5. mutation_gene subset load

- `loadMolecularProfiles(..., feature_type = "mutation_gene", subset, format = "wide")`
  - original: `0.051 s`
  - wrapper: `0.017 s`
  - wrapper warm: `0.013 s`

Interpretation:

- discrete data also benefits from reduced repeated metadata and query overhead
- improvement is meaningful but not as dramatic as wide-table sample subset loading

### 6. treatment response subset load

- `loadTreatmentResponse(..., subset)`
  - original: `0.003 s`
  - wrapper: `0.008 s`
  - wrapper warm: `0.004 s`

Interpretation:

- treatment response is not currently a major performance problem
- there is no strong reason to prioritize optimization here

### 7. Matrix write benchmark

- `storeMatricesInDatabase`
  - original: `0.491 s`
  - wrapper: `0.569 s`

Interpretation:

- this wrapper round did not improve write speed
- current evidence supports prioritizing read-path optimization, not write-path changes

## Main Conclusions

### Conclusion 1

The Linux slowdown is not best explained by `ulimit`.

Evidence:

- the dominant cost is in large SQLite read patterns, especially wide matrix loading
- the strongest gains came from query strategy and SQLite interaction changes, not process limits

### Conclusion 2

The main bottleneck is `loadMolecularProfiles()` on wide omics tables, especially sample-subset reads.

Evidence:

- `CCLE_mRNA` full load: about `3.2x` faster with the optimized wrapper
- `CCLE_mRNA` sample subset load: about `60x` faster with the optimized wrapper

### Conclusion 3

The current read path likely wastes time by reading much more data than needed.

Most likely causes:

- repeated schema introspection
- repeated table listing
- repeated distinct-value scans
- wide-table `SELECT *` patterns
- trimming selected samples only after materializing large data in R

### Conclusion 4

Metadata initialization and treatment response reads are not the highest-priority optimization targets.

## Practical Recommendation

Priority should be:

1. optimize `loadMolecularProfiles()` for wide omics tables
2. specifically optimize sample-subset loading
3. reduce repeated SQLite metadata queries
4. only then consider secondary areas such as write speed

## Suggested Actions for DROMA_Set

If internal package changes are allowed in a later round, the most promising changes are:

### High priority

- reuse connections within a workflow instead of reconnecting repeatedly
- cache `dbListTables()` results per database session
- cache `PRAGMA table_info(table)` results per table
- cache repeated `SELECT DISTINCT` results where safe
- for wide matrix tables, select only `feature_id` plus requested sample columns
- avoid loading the full matrix when only a sample subset is needed

### Medium priority

- audit discrete table access to reduce repeated validation queries
- review whether row-count queries are necessary on hot paths
- profile `CCLE_meth` and other frequently used omics tables in the same way as `CCLE_mRNA`

### Lower priority for now

- treatment response optimization
- write-path optimization
- `ulimit` tuning unless there is direct evidence of OS-level limit failures

## Recommendation About `ulimit`

Current judgment:

- `ulimit` is not the primary fix for this problem
- it is only worth pursuing if there is direct evidence such as:
  - `too many open files`
  - memory limit errors
  - process kill due to explicit resource limits

Without that evidence, performance work should stay focused on SQLite access strategy.

## Final Summary

This benchmark strongly suggests that the observed Linux slowdown is mostly due to current SQLite interaction patterns rather than Linux itself.

The most important finding is:

- sample-subset loading from wide omics tables is dramatically more expensive than it needs to be under the current logic

The external benchmark wrapper demonstrates that this can be improved substantially without changing analysis semantics.
