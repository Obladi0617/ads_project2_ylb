# ADS Project 2: Dijkstra with Four Priority Queues

This project compares Dijkstra's shortest-path algorithm with four min-priority-queue implementations:

- Binary heap
- Ternary heap
- Fibonacci heap
- Binomial heap

The implementation is written in C style in src/1.c and uses only standard C headers, so it builds either with g++ in C++17 mode or with gcc in C11 mode.

## Directory layout

- report_dijkstra_heaps.tex and report_dijkstra_heaps.pdf: English report source and compiled PDF, following the Project 1 layout.
- src/1.c: graph loading, four heaps, Dijkstra, benchmark timing, and heap-operation counters.
- data/dimacs/: the NY, COL, and FLA road-network inputs and their provenance.
- scripts/plot_benchmark.py: plots the measured runtimes.
- results/: benchmark CSV files and the runtime plot.
- build/: the locally compiled executable.
- docs/实验记录.md: project history, current measurements, and handoff notes.
- docs/BENCHMARK_NOTES.md: benchmark methodology and detailed analysis.
- docs/heap_performance_analysis.md: measured comparison of all four heaps, with the supporting timing and operation-count data.

## Build and run

Run these commands from the project root:

    g++ -std=c++17 -O2 .\src\1.c -o .\build\dijkstra_benchmark.exe
    .\build\dijkstra_benchmark.exe

With MinGW gcc, compile as C11 and enable ANSI stdio for the `%lld` format used by the CSV writer:

    gcc -std=c11 -O2 -D__USE_MINGW_ANSI_STDIO=1 .\src\1.c -o .\build\dijkstra_benchmark.exe

The benchmark runs 1,000 source-target queries on each of NY, COL, and FLA with all four heaps. It writes the timing and high-level work CSV files under results/. Query answers are recorded so the implementations can be checked against one another.

To regenerate the runtime plot from the timing CSV:

    python .\scripts\plot_benchmark.py

The plot is saved as results/benchmark_runtime.png.

## Benchmark output files

- results/benchmark_results.csv: total and average runtime for each heap and dataset, plus the number of matching answers.
- results/benchmark_queries.csv: the 1,000 query pairs per dataset and the distance returned by each heap.
- results/benchmark_operations.csv: high-level Dijkstra work, including scanned arcs, settled vertices, inserts, extract-min operations, and successful decrease-key operations.
- results/heap_internal_stats.csv: detailed heap-operation counts, described below.
- results/benchmark_repetitions.csv: raw timings from all three measured repetitions.
- results/benchmark_statistics.csv: median, mean, minimum, maximum, and sample standard deviation for each dataset and heap.
- results/benchmark_environment.json: machine, compiler, build flags, and input/source hashes for the local reproduction.
- results/benchmark_reproduction.log: output from the complete warm-up and three accepted measured rounds.

## Reproduce the repeated experiment on macOS

With Python 3.11 or newer and Apple Clang installed, run:

    python3 scripts/reproduce_benchmark.py
    python3 scripts/plot_benchmark.py

The reproduction script compiles the unchanged C implementation with `-std=c11 -O2 -Wall -Wextra`. It runs one complete warm-up and three sequential measured rounds. Each round uses the same 1,000 distinct query pairs per graph. Warm-up timings are excluded from the summary. The script checks query uniqueness, agreement among all four answers, and identical queries and high-level work across rounds before publishing results. The existing `heap_internal_stats.csv` remains legacy counting data and is explicitly labeled that way in the report; it was not recollected in this local timing experiment.

`benchmark_results.csv` stores the per-heap median of the three measured batch times. The runtime plot uses those medians with observed minimum-to-maximum whiskers; these are not confidence intervals. Temporary binaries and individual pass files remain under the ignored `build/` directory. This full experiment runs sequentially and can take tens of minutes. Completed passes can be revalidated and summarized without rerunning them using `python3 scripts/reproduce_benchmark.py --summarize-existing build/reproduce-<run-id>`.

An independent correctness check is available with `python3 scripts/check_shortest_paths.py`. It compares 3,000 queries on three synthetic 40-vertex graphs against a separate Python Floyd-Warshall implementation, including unreachable pairs, zero-weight/parallel edges, and distances above 32-bit range. These fixtures are separate from the USA road-network timing inputs. The summary is written to `results/reference_checks.json`.

## Counting heap operations

Heap-operation counters are compiled out of the timing build. To collect them, build the separate counting binary and run its count-only pass:

    g++ -std=c++17 -O2 -DCOUNT_HEAP_OPS .\src\1.c -o .\build\dijkstra_benchmark_count.exe
    .\build\dijkstra_benchmark_count.exe --count-heap-ops

The timing binary rejects `--count-heap-ops` and asks for the `-DCOUNT_HEAP_OPS` build.

The counting pass checks the answers again and writes results/heap_internal_stats.csv. It does not overwrite the timing CSV files. Counter collection is separate from the timed pass, so its overhead is not included in benchmark_results.csv.

Each row in heap_internal_stats.csv summarizes one heap on one dataset across the 1,000 queries. The columns mean:

- dataset, heap, queries: dataset name, heap implementation, and number of queries.
- scanned_arcs, settled_vertices: Dijkstra-level work totals.
- pq_inserts, pq_extracts, pq_decrease_keys: priority-queue operation totals.
- comparisons: calls to the key/item comparison function. Binomial root-list degree comparisons are not included.
- upward_steps: parent checks while moving an item upward in an array heap or during Binomial decrease-key.
- downward_iterations: sift-down loop iterations, including the final iteration that stops movement.
- array_swaps: swaps in the binary and ternary heap arrays.
- key_swaps: item/key swaps while Binomial decrease-key moves an item toward the root.
- tree_links: times one tree root is attached beneath another.
- cuts: Fibonacci-heap cuts.
- cascade_steps: parent checks during Fibonacci cascading cuts.
- root_list_steps: implementation-specific traversal work. Fibonacci counts roots copied during consolidation; Binomial counts root-list merge, consolidation, and minimum-root scan iterations. These are useful per implementation but are not identical events.
- children_promoted: children moved to the root list after extract-min.
- max_queue_size: largest priority-queue size observed.
- max_root_count: largest root-list size observed for the linked heaps.
- max_degree: largest number of direct children observed on a node; this is not a general tree-height measurement.

The operation counters help explain why structures with the same Dijkstra-level work can have different runtimes. Compare the detailed counts together with the timings; counters describe work performed and do not measure cache misses or pointer-update latency directly.

## Road-network data

The NY, COL, and FLA DIMACS .gr files are included under data/dimacs so the benchmark can be run after cloning. Their source, reconstruction history, official download page, and graph sizes are documented in data/dimacs/README.md.

The assignment asks the final report or submission package to list dataset download links instead of attaching the raw datasets. This repository keeps the data for development; prepare the course submission according to that instruction.
