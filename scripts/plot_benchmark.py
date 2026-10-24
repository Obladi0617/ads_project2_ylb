import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent.parent
csv_path = root / "results" / "benchmark_results.csv"
image_path = root / "results" / "benchmark_runtime.png"

with csv_path.open("r", newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

if not rows:
    raise SystemExit("benchmark_results.csv contains no results")

arcs_millions = [int(row["arcs"]) / 1_000_000 for row in rows]
binary = [float(row["binary_seconds"]) for row in rows]
ternary = [float(row["ternary_seconds"]) for row in rows]
fibonacci = [float(row["fibonacci_seconds"]) for row in rows]
binomial = [float(row["binomial_seconds"]) for row in rows]

fig, ax = plt.subplots(figsize=(8.5, 5.2))
ax.plot(arcs_millions, binary, marker="o", linewidth=2,
        label="Binary heap")
ax.plot(arcs_millions, ternary, marker="D", linewidth=2,
        label="Ternary heap")
ax.plot(arcs_millions, fibonacci, marker="s", linewidth=2,
        label="Fibonacci heap")
ax.plot(arcs_millions, binomial, marker="^", linewidth=2,
        label="Binomial heap")

ax.set_xticks(arcs_millions,
              [f"{row['dataset']}\n{int(row['arcs']):,}"
               for row in rows])
ax.set_xlabel("Graph size (number of directed arcs)")
ax.set_ylabel("CPU time for 1,000 queries (seconds)")
ax.set_title("Dijkstra runtime by graph size")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(image_path, dpi=300)
print("Saved", image_path)
