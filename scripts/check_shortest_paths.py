"""Check the real benchmark executable against an independent Floyd-Warshall oracle."""

import csv
import hashlib
import json
from pathlib import Path
import random
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parent.parent
N = 40
HEAPS = ("binary", "ternary", "fibonacci", "binomial")


def fixtures():
    rng = random.Random(20261025)
    connected = [(u, (u + 1) % N, rng.randrange(1, 21)) for u in range(N)]
    connected += [((u + 1) % N, u, rng.randrange(1, 21)) for u in range(N)]
    connected += [(rng.randrange(N), rng.randrange(N), rng.randrange(21)) for _ in range(100)]
    disconnected = []
    for start in (0, 20):
        disconnected += [(start + u, start + (u + 1) % 20, 1) for u in range(20)]
        disconnected += [(start + rng.randrange(20), start + rng.randrange(20),
                          rng.randrange(21)) for _ in range(50)]
    disconnected += [(0, 1, 0), (0, 1, 9), (25, 25, 0)]
    large = [(u, u + 1, 1_000_000_000) for u in range(34)]
    large += [(0, 1, 1_500_000_000), (10, 10, 0)]
    return {"NY": connected, "COL": disconnected, "FLA": large}


def floyd_warshall(edges):
    infinity = 10**30
    distance = [[infinity] * N for _ in range(N)]
    for u in range(N):
        distance[u][u] = 0
    for u, v, weight in edges:
        distance[u][v] = min(distance[u][v], weight)
    for k in range(N):
        for u in range(N):
            if distance[u][k] == infinity:
                continue
            for v in range(N):
                distance[u][v] = min(distance[u][v], distance[u][k] + distance[k][v])
    return [[-1 if value == infinity else value for value in row] for row in distance]


def main():
    (ROOT / "build").mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="reference-", dir=ROOT / "build"))
    data = work / "data/dimacs"
    data.mkdir(parents=True)
    (work / "results").mkdir()
    executable = work / "benchmark"
    subprocess.run(["clang", "-std=c11", "-O2", "-Wall", "-Wextra",
                    str(ROOT / "src/1.c"), "-o", str(executable)], check=True)
    oracle = {}
    for name, edges in fixtures().items():
        with (data / f"USA-road-d.{name}.gr").open("w") as stream:
            stream.write(f"p sp {N} {len(edges)}\n")
            for u, v, weight in edges:
                stream.write(f"a {u + 1} {v + 1} {weight}\n")
        oracle[name] = floyd_warshall(edges)
    subprocess.run([str(executable)], cwd=work, check=True, capture_output=True, text=True)
    with (work / "results/benchmark_queries.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    counts = {name: {"queries": 0, "unreachable": 0, "above_int32": 0} for name in oracle}
    pairs = {name: set() for name in oracle}
    for row in rows:
        name = row["dataset"]
        u, v = int(row["source"]) - 1, int(row["target"]) - 1
        if u == v or (u, v) in pairs[name]:
            raise RuntimeError("Benchmark emitted an invalid or duplicate query")
        pairs[name].add((u, v))
        expected = oracle[name][u][v]
        for heap in HEAPS:
            if int(row[f"{heap}_distance"]) != expected:
                raise RuntimeError(f"Incorrect {heap} answer on {name}: {u} -> {v}")
        counts[name]["queries"] += 1
        counts[name]["unreachable"] += expected == -1
        counts[name]["above_int32"] += expected > 2**31 - 1
    if any(group["queries"] != 1000 for group in counts.values()):
        raise RuntimeError("Expected 1,000 query pairs per synthetic fixture")
    if not counts["COL"]["unreachable"] or not counts["FLA"]["above_int32"]:
        raise RuntimeError("Required unreachable/64-bit cases were not exercised")
    result = {"oracle": "Floyd-Warshall, independent Python implementation",
              "vertices_per_fixture": N, "fixture_seed": 20261025,
              "source_sha256": hashlib.sha256((ROOT / "src/1.c").read_bytes()).hexdigest(),
              "query_pairs_checked": len(rows), "heap_answers_checked": len(rows) * 4,
              "input_kind": "synthetic graphs, separate from the USA road-network benchmark",
              "fixtures": {"connected": counts["NY"], "disconnected": counts["COL"],
                           "long_paths_and_isolated_vertices": counts["FLA"]},
              "status": "passed"}
    (ROOT / "results/reference_checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
