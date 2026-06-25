import subprocess, sys

cmds = [
    [sys.executable, "scripts/audit_dataset.py"],
    [sys.executable, "train.py"],
    [sys.executable, "tests/smoke_test_api.py"],
]

for cmd in cmds:
    print("\nRUN:", " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise SystemExit(result.returncode)
