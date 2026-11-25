"""Render the authoritative PDFs so PNG previews have identical text layout."""
from pathlib import Path
import subprocess

figures = Path('figures')
for name in ['recovery-summary', 'neff-through-time', 'p-value-percentages',
             'p-value-distribution', 'trial-001-overview',
             'scenario-p-value-percentages', 'scenario-p-value-distribution',
             'scenario-neff-through-time']:
    subprocess.run(['pdftoppm', '-f', '1', '-singlefile', '-r', '300', '-png',
                    str(figures / (name + '.pdf')), str(figures / name)], check=True)
subprocess.run(['pdftoppm', '-f', '1', '-singlefile', '-r', '300', '-png',
                str(figures / 'trial-001-unit-acceleration.pdf'),
                str(figures / 'trial-001-unit-01')], check=True)
