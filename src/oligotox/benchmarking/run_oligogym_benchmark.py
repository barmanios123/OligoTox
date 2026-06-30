from pathlib import Path
import yaml

def main():
    repo = Path(__file__).resolve().parents[3]
    cfg_path = repo / "benchmarks" / "oligogym" / "configs" / "benchmark_plan.yaml"
    cfg = yaml.safe_load(cfg_path.read_text())
    print("Loaded benchmark plan:", cfg["benchmark_name"])
    print("Datasets planned:", [d["name"] for d in cfg["datasets"]])
    print("Models planned:", [m["name"] for m in cfg["models"]])

if __name__ == "__main__":
    main()
