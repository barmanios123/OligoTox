import importlib

PLATFORMS = [
    # Core registry surfaces
    "oligotox.registry",
    "oligotox.models.registry",
    "oligotox.omics.registry",
    "oligotox.efficacy.registry",
    # Any other notable integration points we can see from the repo layout
    "oligotox.models.efficacy.external_models",
    "oligotox.partner_import",
    "oligotox.evaluation.benchmark_registry",
]

def main():
    print("=== ML platform import check ===")
    for name in PLATFORMS:
        try:
            m = importlib.import_module(name)
            print(f"[OK] {name} imported as {m!r}")
        except Exception as e:
            print(f"[FAIL] {name} -> {e.__class__.__name__}: {e}")

if __name__ == "__main__":
    main()
