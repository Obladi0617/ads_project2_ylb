"""Reproduce three timed rounds plus a complete warm-up.

Run on macOS from any directory with Python 3 and clang installed. Each invocation of
the unchanged benchmark uses 1,000 fixed query pairs per road network. Results
are published only after answers and high-level work match across every pass.
"""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import tempfile
from datetime import datetime, timezone


ROOT = Path(__file__).resolve().parent.parent
HEAPS = ("binary", "ternary", "fibonacci", "binomial")
DATASETS = ("NY", "COL", "FLA")
REPETITIONS = 3


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def command_text(command):
    return subprocess.check_output(command, text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summarize-existing", type=Path,
                        help="validate and summarize completed passes without rerunning them")
    args = parser.parse_args()
    build = ROOT / "build"
    build.mkdir(exist_ok=True)
    work = (args.summarize_existing.resolve() if args.summarize_existing else
            Path(tempfile.mkdtemp(prefix="reproduce-", dir=build)))
    timed = work / "benchmark"
    flags = ["-std=c11", "-O2", "-Wall", "-Wextra"]
    if not args.summarize_existing:
        subprocess.run(["clang", *flags, str(ROOT / "src/1.c"), "-o", str(timed)], check=True)
    environment = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "macos_version": platform.mac_ver()[0],
        "architecture": platform.machine(),
        "processor": command_text(["sysctl", "-n", "machdep.cpu.brand_string"]),
        "memory_bytes": int(command_text(["sysctl", "-n", "hw.memsize"])),
        "physical_cores": int(command_text(["sysctl", "-n", "hw.physicalcpu"])),
        "logical_cores": os.cpu_count(),
        "compiler": command_text(["clang", "--version"]),
        "build_flags": flags,
        "source_commit": command_text(["git", "-C", str(ROOT), "rev-parse", "HEAD"]),
        "source_sha256": digest(ROOT / "src/1.c"),
        "timed_binary_sha256": digest(timed),
        "warmup_rounds": 1,
        "measured_rounds": REPETITIONS,
        "queries_per_dataset_per_round": 1000,
        "clock": "clock(): process CPU seconds; loading and allocation excluded",
        "heap_order": "dataset-dependent rotation from src/1.c; same order on each repetition",
        "internal_counters": "pre-existing legacy CSV retained; not remeasured in this experiment",
        "excluded_partial_rounds": [4] if (work / "round_04").exists() else [],
        "input_sha256": {name: digest(ROOT / f"data/dimacs/USA-road-d.{name}.gr")
                         for name in DATASETS},
    }
    manifest = work / "benchmark_environment.json"
    if args.summarize_existing:
        recorded = json.loads(manifest.read_text())
        for field in ("source_sha256", "timed_binary_sha256", "input_sha256"):
            if recorded[field] != environment[field]:
                raise RuntimeError(f"Recorded experiment no longer matches {field}")
        environment = recorded
        environment["measured_rounds"] = REPETITIONS
        environment["excluded_partial_rounds"] = [4] if (work / "round_04").exists() else []
    else:
        manifest.write_text(json.dumps(environment, indent=2) + "\n", encoding="utf-8")
    reference = None
    measurements = []
    log_path = work / "benchmark_reproduction.log"

    def run_pass(label, binary, arguments=()):
        folder = work / label
        folder.mkdir()
        (folder / "data").symlink_to(ROOT / "data", target_is_directory=True)
        (folder / "results").mkdir()
        marker = f"{datetime.now(timezone.utc).isoformat()} START {label}"
        print(marker, flush=True)
        with log_path.open("a", encoding="utf-8") as log:
            log.write(marker + "\n")
            with subprocess.Popen([str(binary), *arguments], cwd=folder, text=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as process:
                for line in process.stdout:
                    print(f"[{label}] {line}", end="", flush=True)
                    log.write(line)
                    log.flush()
                return_code = process.wait()
            if return_code:
                raise RuntimeError(f"{label} failed with exit code {return_code}; see {log_path}")
            log.write(f"END {label}\n")
        return folder / "results"

    for repetition in range(REPETITIONS + 1):
        label = "warmup" if repetition == 0 else f"round_{repetition:02d}"
        output = (work / label / "results" if args.summarize_existing else
                  run_pass(label, timed))
        rows = read_csv(output / "benchmark_results.csv")
        if tuple(row["dataset"] for row in rows) != DATASETS:
            raise RuntimeError(f"Unexpected dataset set in {label}")
        if any(int(row["queries"]) != 1000 or int(row["matching_answers"]) != 1000
               for row in rows):
            raise RuntimeError(f"Incomplete or mismatching queries in {label}")
        queries = read_csv(output / "benchmark_queries.csv")
        for name in DATASETS:
            group = [row for row in queries if row["dataset"] == name]
            pairs = {(row["source"], row["target"]) for row in group}
            if len(group) != 1000 or len(pairs) != 1000:
                raise RuntimeError(f"Query pairs are incomplete or repeated in {label}/{name}")
            for row in group:
                if row["source"] == row["target"] or len({row[f"{h}_distance"] for h in HEAPS}) != 1:
                    raise RuntimeError(f"Invalid query or unequal answers in {label}/{name}")
        hashes = {name: digest(output / name) for name in
                  ("benchmark_queries.csv", "benchmark_operations.csv")}
        if reference is None:
            reference = hashes
        elif hashes != reference:
            raise RuntimeError(f"Queries, answers, or high-level work changed in {label}")
        if repetition:
            measurements.extend({"repetition": repetition, **row} for row in rows)

    publish = work / "published"
    publish.mkdir(exist_ok=True)
    fields = list(measurements[0])
    write_csv(publish / "benchmark_repetitions.csv", fields, measurements)
    summaries = []
    statistics_rows = []
    for name in DATASETS:
        group = [row for row in measurements if row["dataset"] == name]
        summary = {key: value for key, value in group[0].items() if key != "repetition"}
        for heap in HEAPS:
            values = [float(row[f"{heap}_seconds"]) for row in group]
            median = statistics.median(values)
            summary[f"{heap}_seconds"] = f"{median:.6f}"
            summary[f"{heap}_avg_ms"] = f"{median:.3f}"
            statistics_rows.append({
                "dataset": name, "heap": heap, "repetitions": REPETITIONS,
                "median_seconds": f"{median:.6f}",
                "mean_seconds": f"{statistics.mean(values):.6f}",
                "min_seconds": f"{min(values):.6f}",
                "max_seconds": f"{max(values):.6f}",
                "sample_stddev_seconds": f"{statistics.stdev(values):.6f}",
            })
        summaries.append(summary)
    write_csv(publish / "benchmark_results.csv", list(summaries[0]), summaries)
    write_csv(publish / "benchmark_statistics.csv", list(statistics_rows[0]), statistics_rows)
    for name in ("benchmark_queries.csv", "benchmark_operations.csv"):
        shutil.copyfile(work / "round_01/results" / name, publish / name)
    environment.update(reference)
    log = log_path.read_text()
    environment["started_utc"] = log.splitlines()[0].split(" START ")[0]
    environment["summarized_utc"] = datetime.now(timezone.utc).isoformat()
    (publish / "benchmark_environment.json").write_text(
        json.dumps(environment, indent=2) + "\n", encoding="utf-8")
    last_marker = f"END round_{REPETITIONS:02d}\n"
    if last_marker not in log:
        raise RuntimeError("Final measured pass did not finish")
    (publish / log_path.name).write_text(log[:log.index(last_marker) + len(last_marker)])
    for path in publish.iterdir():
        shutil.copyfile(path, ROOT / "results" / path.name)
    print(f"Published {REPETITIONS}-round medians and validated data; workspace: {work}", flush=True)


if __name__ == "__main__":
    main()
