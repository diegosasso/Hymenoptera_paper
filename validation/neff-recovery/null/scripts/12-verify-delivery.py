"""Hash portable cluster outputs, or verify their transfer and preserved inputs."""
from pathlib import Path
import csv
import hashlib
import json
import sys

root = Path.cwd()

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

if '--hash-only' in sys.argv:
    hashes = {}
    for directory in ['simulated', 'results']:
        for trial in range(1, 101):
            folder = root / directory / f'trial-{trial:03d}'
            for path in sorted(folder.iterdir()):
                if path.is_file() and not path.name.endswith('.tmp'):
                    hashes[str(path.relative_to(root))] = digest(path)
    (root / 'provenance/portable-output-sha256.json').write_text(
        json.dumps(hashes, indent=2) + '\n')
    print(f'Hashed {len(hashes)} portable scientific files.')
    sys.exit(0)

messages = []
hashes = json.loads((root / 'provenance/portable-output-sha256.json').read_text())
for name, expected in hashes.items():
    assert digest(root / name) == expected, name
messages.append(f'All {len(hashes)} portable files match Roihu SHA-256 hashes.')

preserved = (root / 'provenance/associated-inputs-sha256.txt').read_text().splitlines()
for line in preserved:
    expected, name = line.split('  ', 1)
    assert digest(root.parent / name) == expected, name
messages.append(f'All {len(preserved)} retained associated files match their original hashes.')

for line in (root / 'provenance/analysis-sha256.txt').read_text().splitlines():
    expected, name = line.split('  ', 1)
    assert digest(root / name) == expected, name
messages.append('Prespecified simulation, inference, metric, audit and seed files are unchanged.')

with (root / 'design/task-manifest.csv').open() as stream:
    manifest = list(csv.DictReader(stream))
assert [int(row['trial']) for row in manifest] == list(range(1, 101))
for row in manifest:
    path = root / 'simulated' / f"trial-{int(row['trial']):03d}" / 'tree.phy'
    assert digest(path) == row['source_tree_sha256']
messages.append('All 100 control tree files match their paired source hashes.')

for table in ['inference-status', 'trials']:
    with (root / 'summary' / (table + '.csv')).open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 100 and len({row['trial'] for row in rows}) == 100
    assert all(row['status'] in ['ok', 'ok_constants_omitted'] for row in rows)
messages.append('All 100 planned trials completed inference and summaries.')

for name in ['scenario-p-value-percentages.pdf', 'scenario-p-value-distribution.pdf',
             'scenario-neff-through-time.pdf', 'all-trial-diagnostics.pdf',
             'trial-001-overview.pdf', 'trial-001-unit-acceleration.pdf']:
    assert (root / 'figures' / name).stat().st_size > 1000
assert (root / 'manuscript/methods-preview.pdf').stat().st_size > 1000
messages.append('Comparison figures, complete trial/unit atlases and methods PDF are present.')
text = '\n'.join(messages) + '\n'
(root / 'provenance/delivery-verification.txt').write_text(text)
print(text)
