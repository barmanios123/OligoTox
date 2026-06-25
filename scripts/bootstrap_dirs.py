from pathlib import Path

dirs = [
    'config', 'data/raw', 'data/interim', 'data/processed',
    'artifacts/models', 'artifacts/reports', 'artifacts/logs',
    'demos', 'tools', 'docs'
]
for d in dirs:
    Path(d).mkdir(parents=True, exist_ok=True)
    print('ok', d)
