"""Verify the complete local delivery against the independent Roihu audit."""
import csv
import hashlib
from pathlib import Path

root = Path(__file__).resolve().parent.parent
with (root / "summary/completion-audit.csv").open() as stream:
    audit = list(csv.DictReader(stream))
assert len(audit) == 1500
assert {int(row["task"]) for row in audit} == set(range(1, 1501))
conditions = {}
checked = 0
for row in audit:
    condition = int(row["condition"])
    conditions[condition] = conditions.get(condition, 0) + 1
    task_dir = root / "results" / f"task-{int(row['task']):04d}"
    for filename, column in [
        ("result.rds", "result_md5"),
        ("maps-full.rds", "full_maps_md5"),
        ("maps-study.rds", "study_maps_md5"),
    ]:
        path = task_dir / filename
        digest = hashlib.md5()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        assert digest.hexdigest() == row[column], f"Transfer mismatch: {path}"
        checked += 1
assert conditions == {i: 100 for i in range(1, 16)}
status = (root / "summary/collection-status.txt").read_text()
assert "Completed: 1500 of 1500 paired trials." in status
report = f"PASS: {checked} local result/map files match the complete Roihu audit.\n"
(root / "provenance/transfer-verification.txt").write_text(report)
print(report, end="")

