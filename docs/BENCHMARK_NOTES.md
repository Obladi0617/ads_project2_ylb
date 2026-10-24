# Dijkstra benchmark: three-run local reproduction

Date: 2026-10-25 (Asia/Shanghai).

## Environment and method

- Apple M1 Max, arm64, 10 CPU cores, 32 GiB RAM; macOS 26.6.2.
- Apple Clang 21.0.0; unchanged C source compiled with -std=c11 -O2 -Wall -Wextra.
- One full warm-up and three complete sequential measured program runs.
- NY, COL, and FLA each use the same 1,000 distinct ordered source-target pairs per run.
- Four heaps per dataset; 36,000 timed query executions. Repetition does not create new distinct pairs.
- Each heap/dataset CPU batch time is measured with clock(); graph loading, query generation, and allocation are excluded. Per-query cleanup is included.
- The starting heap rotates by dataset, while its order stays fixed across repetitions.
- benchmark_results.csv now stores medians of the three complete runs. The former single-run Windows times are historical Git data, not inputs to these statistics.
- A fourth run had started before the user requested three runs. It was terminated; none of its values enters the summary. The warm-up is also excluded.

## Median results

| Dataset | Binary (s) | Ternary (s) | Fibonacci (s) | Binomial (s) |
| --- | --- | --- | --- | --- |
| NY | 12.744 | 11.608 | 24.264 | 21.850 |
| COL | 20.584 | 18.823 | 37.767 | 33.708 |
| FLA | 51.513 | 47.112 | 96.150 | 85.851 |

Ternary was faster than binary in each of the three paired observations on each network. All four heaps matched on 1,000/1,000 pairs per graph in every measured run. Query/answer files and high-level Dijkstra operation files are byte-identical across the warm-up and the three measured runs.

## Observed variability

| Dataset | Heap | Median (s) | Min (s) | Max (s) | Sample SD (s) |
| --- | --- | --- | --- | --- | --- |
| NY | binary | 12.744330 | 12.687818 | 12.766622 | 0.040621 |
| NY | ternary | 11.607607 | 11.569196 | 11.646938 | 0.038872 |
| NY | fibonacci | 24.264308 | 24.263044 | 24.266950 | 0.001993 |
| NY | binomial | 21.850092 | 21.822269 | 21.904328 | 0.041732 |
| COL | binary | 20.584066 | 20.481893 | 20.661779 | 0.090220 |
| COL | ternary | 18.822622 | 18.786543 | 18.949682 | 0.085694 |
| COL | fibonacci | 37.767477 | 37.620855 | 37.967750 | 0.174138 |
| COL | binomial | 33.708046 | 33.631397 | 33.759263 | 0.064353 |
| FLA | binary | 51.513346 | 51.369597 | 51.790087 | 0.213721 |
| FLA | ternary | 47.111732 | 47.066759 | 47.962900 | 0.504906 |
| FLA | fibonacci | 96.149886 | 95.915817 | 96.384774 | 0.234479 |
| FLA | binomial | 85.850519 | 85.569537 | 86.079252 | 0.255303 |

The runtime plot uses medians with observed minimum-to-maximum whiskers. These are descriptive ranges, not confidence intervals.

## Independent correctness check

scripts/check_shortest_paths.py compares the actual executable against an independent Python Floyd-Warshall oracle on three synthetic 40-vertex graphs, separate from the USA inputs. All 12,000 heap answers for 3,000 pairs passed, including 1,122 unreachable pairs and 330 distances beyond signed 32-bit range. Cases also include zero weights, parallel edges, and isolated vertices. The output is results/reference_checks.json.

## Interpretation and limits

The high-level graph work is identical for all heaps. Three repetitions describe this fixed workload on this machine; they do not establish a universal ranking or sample different query distributions. CPU frequency, thermal state, and background load were not locked.

The retained heap_internal_stats.csv comes from the earlier implementation's independent counter pass. It was not regenerated during the three-run local experiment. Use it as legacy evidence of structural patterns, not as an exact decomposition of the new CPU times. In particular, its root-list step definitions differ for Fibonacci and binomial, and binomial degree comparisons are excluded from the key/item comparison field.

Local graphs were reconstructed from public edge arrays and checked against published vertex/arc counts. They have not been compared arc by arc with official downloads; the official HTTPS host still timed out during this work. The report lists official download links and does not include raw networks in its submission materials.

## Reproduction and records

Run python3 scripts/reproduce_benchmark.py, then python3 scripts/plot_benchmark.py. Reproduce the oracle with python3 scripts/check_shortest_paths.py. Completed pass directories can be revalidated with --summarize-existing.

Raw values: benchmark_repetitions.csv. Summary: benchmark_results.csv and benchmark_statistics.csv. Exact inputs and compiler/source hashes: benchmark_environment.json. Complete accepted-pass output: benchmark_reproduction.log. Query answers and high-level work: benchmark_queries.csv and benchmark_operations.csv.
