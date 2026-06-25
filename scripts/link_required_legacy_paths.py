from pathlib import Path

ROOT = Path.cwd()

ROOT.joinpath("archive/__init__.py").touch()
ROOT.joinpath("archive/app_mvp").mkdir(parents=True, exist_ok=True)
ROOT.joinpath("archive/app_mvp/__init__.py").touch()

ROOT.joinpath("app_mvp").mkdir(parents=True, exist_ok=True)
ROOT.joinpath("app_mvp/__init__.py").touch()

mvp = ROOT / "app_mvp/mvp_support.py"
if not mvp.exists():
    mvp.write_text("from archive.app_mvp.mvp_support import *\n")

for name in ["omics_study_manifest_demo.csv", "omics_study_manifest_template.csv"]:
    legacy = ROOT / "archive/_legacy/manifests" / name
    live = ROOT / "manifests" / name
    live.parent.mkdir(parents=True, exist_ok=True)
    if legacy.exists() and not live.exists():
        live.write_bytes(legacy.read_bytes())

print("Legacy path linking complete.")
