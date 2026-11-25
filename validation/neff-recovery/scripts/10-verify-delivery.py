"""Verify the complete counts and the values used in the delivered tables."""
from pathlib import Path
import csv
import hashlib
import json
import math

root = Path('.')
def rows(name):
    with (root / name).open() as stream:
        return list(csv.DictReader(stream))

trials = rows('summary/trials.csv')
correlations = rows('summary/correlations.csv')
percentages = rows('summary/p-value-percentages.csv')
inference = rows('summary/inference-status.csv')
audit = rows('provenance/metric-audit.csv')
assert len(trials) == len(inference) == len(audit) == 100
assert {int(r['trial']) for r in trials} == set(range(1, 101))
assert all(r['status'] in ('ok', 'ok_constants_omitted') for r in inference)
assert all(int(r['n_fit_failed']) == int(r['n_map_failed']) == 0 for r in inference)
assert sum(int(r['n_maps']) for r in inference) == sum(int(r['n_variable']) for r in trials) * 100
assert all(r['status'] == 'passed' and float(r['maximum_absolute_error']) < 1e-10 for r in audit)
for result in percentages:
    matching = [r for r in correlations if r['source'] == result['source'] and r['cutoff'] == result['cutoff']]
    valid = [r for r in matching if r['nominal_p'] != 'NA' and r['Pearson_r'] != 'NA']
    detected = [r for r in valid if float(r['nominal_p']) < float(result['alpha'])]
    if result['direction'] == 'positive_r':
        detected = [r for r in detected if float(r['Pearson_r']) > 0]
    assert int(result['n_estimable']) == len(valid)
    assert int(result['n_detected']) == len(detected)
    assert float(result['percent_all_trials']) == len(detected)
    assert int(result['n_unavailable']) == 100 - len(valid)
for line in Path('provenance/analysis-sha256.txt').read_text().splitlines():
    digest, name = line.split('  ', 1)
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == digest, name
for item in json.loads(Path('provenance/ontophylo-source-sha256.json').read_text()):
    file = Path('provenance/ontophylo') / item['file']
    assert hashlib.sha256(file.read_bytes()).hexdigest() == item['sha256']
for name in ['main-text-methods-standalone.tex', 'supplementary-methods-standalone.tex']:
    text = (Path('manuscript') / name).read_text()
    assert '\\Neff' not in text, 'Unexpanded result macro in '+name
required = ['recovery-summary.pdf', 'neff-through-time.pdf', 'p-value-percentages.pdf',
            'p-value-distribution.pdf', 'all-trial-neff-curves.pdf',
            'all-trial-diagnostics.pdf', 'trial-001-overview.pdf',
            'trial-001-unit-acceleration.pdf']
assert all((Path('figures') / name).is_file() for name in required)
assert Path('manuscript/methods-preview.pdf').is_file()
report = [f'Verified {len(trials)} trials and {sum(int(r["n_maps"]) for r in inference)} maps.',
          'All variable-character fits and maps succeeded.',
          f'Independent Neff checks: {sum(int(r["independent_neff_values"]) for r in audit)} values.',
          'P-value counts independently reproduced from per-trial CSVs.',
          'Frozen inference code, seed manifest and 91 package files unchanged.',
          'Literal-number LaTeX snippets, compiled document and required figures present.']
Path('provenance/delivery-verification.txt').write_text('\n'.join(report)+'\n')
print('\n'.join(report))
