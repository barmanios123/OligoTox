from __future__ import annotations

import os
import sys
from pathlib import Path


def bootstrap_project() -> Path:
    cwd = Path.cwd().resolve()

    if cwd.name == "downloads":
        project_root = cwd.parent
    elif (cwd / "oligotox_app.py").exists():
        project_root = cwd
    else:
        probe = cwd
        while probe != probe.parent:
            if (probe / "oligotox_app.py").exists():
                project_root = probe
                break
            probe = probe.parent
        else:
            project_root = cwd

    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    os.chdir(project_root)
    return project_root


if __name__ == "__main__":
    root = bootstrap_project()
    print("Bootstrapped project root:", root)
