import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent.parent
csv_path = root / "results" / "benchmark_results.csv"
image_path = root / "results" / "benchmark_runtime.png"
statistics_path = root / "results" / "benchmark_statistics.csv"

with csv_path.open("r", newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

if not rows:
    raise SystemExit("benchmark_results.csv contains no results")

arcs_millions = [int(row["arcs"]) / 1_000_000 for row in rows]
binary = [float(row["binary_seconds"]) for row in rows]
ternary = [float(row["ternary_seconds"]) for row in rows]
fibonacci = [float(row["fibonacci_seconds"]) for row in rows]
binomial = [float(row["binomial_seconds"]) for row in rows]

with statistics_path.open("r", newline="", encoding="utf-8") as f:
    statistics_rows = {(row["dataset"], row["heap"]): row
                       for row in csv.DictReader(f)}


def observed_range(heap, medians):
    """Asymmetric whiskers show observed min/max, not confidence intervals."""
    lower, upper = [], []
    for row, median in zip(rows, medians):
        stats = statistics_rows[(row["dataset"], heap)]
        lower.append(median - float(stats["min_seconds"]))
        upper.append(float(stats["max_seconds"]) - median)
    return lower, upper

fig, ax = plt.subplots(figsize=(8.5, 5.2))
for heap, values, marker, label in (
        ("binary", binary, "o", "Binary heap"),
        ("ternary", ternary, "D", "Ternary heap"),
        ("fibonacci", fibonacci, "s", "Fibonacci heap"),
        ("binomial", binomial, "^", "Binomial heap")):
    ax.errorbar(arcs_millions, values, yerr=observed_range(heap, values),
                marker=marker, linewidth=2, capsize=4, label=label)

ax.set_xticks(arcs_millions,
              [f"{row['dataset']}\n{int(row['arcs']):,}"
               for row in rows])
ax.set_xlabel("Graph size (number of directed arcs)")
ax.set_ylabel("Median CPU time for 1,000 queries (seconds)")
ax.set_title("Dijkstra runtime: three repetitions; whiskers show min-max")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(image_path, dpi=300)
print("Saved", image_path)
