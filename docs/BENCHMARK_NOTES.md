# Dijkstra heap benchmark notes

## Results

Each dataset uses 1,000 deterministic point-to-point queries. All four heaps returned identical shortest-path distances for all queries.

| Dataset | Vertices | Arcs | Binary (s) | Ternary (s) | Fibonacci (s) | Binomial (s) | Matching |
|---|---:|---:|---:|---:|---:|---:|---:|
| NY | 264,346 | 733,846 | 9.752 | 9.929 | 28.822 | 22.326 | 1000/1000 |
| COL | 435,666 | 1,057,066 | 15.880 | 16.062 | 45.283 | 34.043 | 1000/1000 |
| FLA | 1,070,376 | 2,712,798 | 40.315 | 41.061 | 115.971 | 87.850 | 1000/1000 |

Average time per query in milliseconds is numerically the same as the total seconds above because each batch contains 1,000 queries.

The ternary heap took 1.1% to 1.9% longer than the binary heap in this run. Fibonacci took 2.85 to 2.96 times as long as binary; binomial took 2.14 to 2.29 times as long. The binary and ternary results are close enough that one run cannot establish which is reliably faster.

## Work performed

results/benchmark_operations.csv records the Dijkstra work counters. Each heap performed the same work on a given network; the deterministic (distance, vertex) tie break kept the extracted-vertex sequence consistent.

| Dataset | Scanned arcs | Settled vertices / extracts | Inserts | Decrease-key |
|---|---:|---:|---:|---:|
| NY | 385,059,749 | 135,502,866 | 135,947,390 | 19,297,309 |
| COL | 541,176,229 | 220,809,144 | 221,168,701 | 16,421,882 |
| FLA | 1,347,518,609 | 531,602,948 | 532,034,751 | 50,489,498 |

These are totals over the 1,000 queries, not counts for one query. The search stops when the target is extracted, so it does not scan the entire graph for each pair. In these workloads, extract-min calls outnumber successful decrease-key calls by about 7 to 13 times. That makes the cost of repeated extraction and its heap maintenance especially relevant; Fibonacci's amortized constant-time decrease-key does not make the whole search faster by itself.

## What the timings suggest

- Binary and ternary heaps both store nodes contiguously in arrays and use a position array for decrease-key. They have good memory locality and avoid per-node pointers.
- A ternary heap has fewer levels: its height is approximately log_3(n) instead of log_2(n), around 63% as many levels. Insert and decrease-key therefore travel through fewer parents.
- On extract-min, a ternary heap may inspect three children at each level, while a binary heap inspects at most two. The shorter path does not automatically mean less work. In this implementation, the extra child comparison and index bookkeeping nearly offset the reduced height; measured times were within 2%.
- A binomial heap maintains linked trees and its root list. Extract-min scans roots, reverses/merges child roots, and links trees of equal degree. Decrease-key can move a node up its parent chain. These pointer updates and traversals are more costly than array index operations here.
- A Fibonacci heap makes insert and decrease-key cheap in its amortized analysis. Extract-min is more involved: it promotes the minimum root's children and consolidates roots of equal degree. Its linked root and child lists require pointer chasing and list updates. Road-query workloads still perform hundreds of millions of extracts, so those costs outweigh the theoretical savings from decrease-key in this implementation.

This comparison measures full Dijkstra query batches, not isolated heap operations. The operation counts show that all heaps processed the same number of graph edges and priority-queue operations; the remaining time difference comes from the cost of carrying out those operations and their memory-access patterns. A hardware profiler or a dedicated heap microbenchmark would be needed to split cache misses, comparisons, and pointer updates precisely.

The binary/ternary gap is only 1% to 2%, and this is one timed run per heap and graph. Treat that gap as inconclusive until repeated runs are compared by their median. The larger gaps between array heaps and linked heaps are consistent across all three graphs, but repeated measurements would also make those estimates stronger.

## Method

- The program loads each DIMACS graph before timing and generates 1,000 unique source-target pairs from the full vertex range; source and target differ.
- All four heaps receive the same query list. results/benchmark_queries.csv stores the pairs and four answers.
- Each point-to-point Dijkstra run stops when the target is extracted as the current minimum.
- clock() measures CPU time around each 1,000-query batch. Graph loading, query generation, and heap/workspace allocation are excluded. Per-query heap and state cleanup is included.
- results/benchmark_operations.csv records high-level Dijkstra counters during the timed runs; all counts match per graph. results/heap_internal_stats.csv comes from a separate count-only pass, so its internal-counter overhead is excluded from the recorded runtimes.
- The benchmark rotates the starting heap across datasets to reduce systematic timing-order bias.
- results/benchmark_results.csv contains total and average query times. results/benchmark_operations.csv contains scanned arcs, settled vertices, and priority-queue operation counts. plot_benchmark.py generates results/benchmark_runtime.png, plotting total time against arc count.

## Internal operation counters

The instrumented pass is run separately with build/dijkstra_benchmark.exe --count-heap-ops. It uses the same query generator and 1,000 pairs per graph, checks the answers, and writes results/heap_internal_stats.csv. It does not replace the timing CSV or include counter collection in the recorded runtime.

The high-level Dijkstra operation counts still match across all heaps. The internal counts show where each structure does different work. NY is a representative example:

| Heap | Key/item comparisons | Upward steps | Downward iterations | Array swaps | Key swaps | Tree links | Cuts | Root/list steps | Children promoted | Max queue | Max roots | Max degree |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Binary | 2,278,648,978 | 192,152,111 | 1,141,483,900 | 1,042,912,603 | 0 | 0 | 0 | 0 | 0 | 1,002 | 0 | 0 |
| Ternary | 2,291,155,185 | 175,210,738 | 822,711,177 | 707,198,507 | 0 | 0 | 0 | 0 | 0 | 1,002 | 0 | 0 |
| Fibonacci | 2,745,491,161 | 0 | 0 | 0 | 0 | 1,020,394,879 | 4,047,134 | 1,686,562,452 | 1,015,908,078 | 1,002 | 22 | 9 |
| Binomial | 1,591,717,381 | 23,300,320 | 0 | 0 | 4,626,161 | 1,037,904,625 | 0 | 4,092,363,948 | 1,037,464,940 | 1,002 | 17 | 9 |

For the array heaps, the ternary heap had 28% fewer downward-loop iterations and 32% fewer swaps than the binary heap. Its upward steps were also lower. Yet its comparator calls were 0.5% higher on NY, 0.9% higher on COL, and 0.8% higher on FLA. This matches the structural tradeoff: fewer levels, with up to three children checked at each downward level instead of two. The measured runtime difference was small, so the timing pass alone does not establish a reliable winner.

The linked heaps performed about 7 to 8 tree links per extract-min across all three graphs. Fibonacci also promoted about 7 to 7.5 children per extract-min. Thus tree linking and child promotion are frequent, substantial parts of their work. The link counts for Fibonacci and binomial are similar, though, while Fibonacci took about 29% longer than binomial on NY. Link count alone does not explain their timing gap; comparisons, list traversal, and pointer updates also matter.

The queue was much smaller than each full graph: its maximum size was 1,002 for NY, 1,037 for COL, and 1,575 for FLA. This reflects point-to-point early termination. The maximum Fibonacci root count was 21 to 23, and the maximum degree was 9 to 10. The binomial root count peaked at 17 to 19, with degree 9 to 10. For a binomial tree, degree gives its height; Fibonacci maximum degree does not give the height of its tallest tree.

Counter definitions and limits:

- comparisons counts key/item comparator calls. Binomial degree comparisons in root-list operations are not included in this column.
- upward_steps counts parent checks for array sift-up and binomial decrease-key. downward_iterations counts sift-down loop iterations, including the final iteration that stops movement.
- array_swaps counts binary/ternary heap swaps; key_swaps counts the vertex/key swaps used by binomial decrease-key.
- tree_links counts calls that attach one tree root below another. cuts and cascade_steps count Fibonacci cuts and parent checks during cascading cuts.
- root_list_steps is implementation-specific: Fibonacci counts roots copied during consolidation; binomial counts iterations in root-list merging, consolidation, and minimum-root scanning. It is useful for understanding each implementation's work, but is not a perfectly identical event across the two heaps.
- children_promoted counts children moved back to the root list after extract-min. max_degree is the largest number of direct children observed, not a general tree-height measurement.

These counters support the conclusion that the linked heaps do substantial structural maintenance, while also showing that the merge/link count by itself is not enough to explain total runtime. A direct measurement of cache misses and pointer-update costs would require hardware profiling. The full per-network, per-heap values are in results/heap_internal_stats.csv.

## Implementations

- Binary heap: array-backed, with a position array for decrease-key.
- Ternary heap: array-backed 3-ary heap, with a position array for decrease-key.
- Fibonacci heap: circular linked root and child lists, consolidation, cuts, and cascading cuts.
- Binomial heap: root list ordered by degree, union/link operations, minimum-root extraction, and decrease-key by moving the item up its parent chain.

## Environment

- Date: 2026-10-23
- CPU: Intel Core i7-14700HX, 20 cores / 28 logical processors
- Compiler: MinGW g++ 15.2.0
- Build flags: g++ -std=c++17 -O2
- The program is single-threaded.

## Data note

The tested local NY, COL, and FLA files are listed in data/dimacs/README.md. That README records that the local graph files were reconstructed from a public repository's graph arrays because the original DIMACS host was unreachable in the earlier download attempt. Mention this provenance in the report.
