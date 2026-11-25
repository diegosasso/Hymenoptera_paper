"""Run independent R fits in a bounded Roihu allocation; preserve each task log."""
import concurrent.futures
import csv
import os
from pathlib import Path
import subprocess
import time

workers = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
Path("logs/tasks").mkdir(parents=True, exist_ok=True)
Path("results").mkdir(exist_ok=True)


def run_task(item):
    stage, task_id = item
    script = "scripts/04-fit-bm.R" if stage == "fit" else "scripts/06-null-comparison.R"
    log_path = Path("logs/tasks") / f"{stage}-{task_id:04d}.log"
    started = time.time()
    with log_path.open("w") as stream:
        result = subprocess.run(
            ["Rscript", "--vanilla", script, str(task_id)],
            stdout=stream, stderr=subprocess.STDOUT, check=False,
        )
    return {"stage": stage, "task_id": task_id, "returncode": result.returncode,
            "elapsed_seconds": round(time.time() - started, 3), "log": str(log_path)}


with open("data/fit-manifest.csv", newline="") as stream:
    fit_ids = [int(row["task_id"]) for row in csv.DictReader(stream)]
stages = [("fit", fit_ids), ("null", list(range(1, 200)))]
failed = 0
with open("results/batch-status.csv", "w", newline="") as status_stream:
    fields = ["stage", "task_id", "returncode", "elapsed_seconds", "log"]
    writer = csv.DictWriter(status_stream, fieldnames=fields)
    writer.writeheader()
    status_stream.flush()
    for stage, task_ids in stages:
        print(f"Starting {stage}: {len(task_ids)} tasks, {workers} workers", flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(run_task, (stage, task_id)) for task_id in task_ids]
            for future in concurrent.futures.as_completed(futures):
                row = future.result()
                writer.writerow(row)
                status_stream.flush()
                failed += row["returncode"] != 0
                print(row, flush=True)
print(f"Finished; failed tasks: {failed}", flush=True)
raise SystemExit(1 if failed else 0)
