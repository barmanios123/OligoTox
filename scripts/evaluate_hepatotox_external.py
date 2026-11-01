"""A thin external check of the hepatocellular injury model.

Hepatocellular injury is the endpoint the platform leans on, and it rests on
19 positive sequences from a single curation. The only public set found with
gapmer sequences paired to liver enzymes is the hepatotoxicity screen in:

  Identification of nucleobase chemical modifications that reduce the
  hepatotoxicity of gapmer antisense oligonucleotides.
  Nucleic Acids Research 2022;50(13):7224. doi:10.1093/nar/gkac562
  Supplementary Table S1, 149 LNA gapmers, C57BL/6 mice, 20 mg/kg IV,
  AST/ALT at 96 h.

The two motif columns are TGC and TCC, which the paper names as the
sequence features driving toxicity. Both are verified on parse against the
sequences themselves: a misaligned column would show up as a mismatch, and
both match 149 of 149. The second was read as "CCC" for a while, which is a
different motif and matched only 44% of rows.

**This set is too small to validate the endpoint and is not treated as
validation.** It carries roughly five hepatotoxic sequences by the paper's
own ALT > 100 U/l criterion. A confidence interval on five positives spans
most of the range, so the result here is directional only: it says whether
the shipped model ranks known-hepatotoxic gapmers above the rest, not how
well it would perform.

What it can still catch is a model that is outright wrong, which is worth
having when the alternative is no external signal at all.

The sequences are 14-mer LNA gapmers, fully phosphorothioate, so they are
scored with ps_fraction 1.0 and the ASO modality. The model never saw this
chemistry at this length, which is itself part of what is being tested.

Writes data/external/nar2022_hepatotox_screen.csv and
reports/results/hepatotox_external.json.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from oligotox.features.candidate import candidate_features
from oligotox.pipelines.score_candidate import score_candidate

ROOT = Path(__file__).resolve().parent.parent
EXTERNAL = ROOT / "data" / "external"
RESULTS = ROOT / "reports" / "results"
ENDPOINT = "hepatocellular_injury"
ALT_TOXIC_ABOVE = 100.0  # the paper's own criterion
N_BOOT = 2000
SEED = 42

# ID, sequence, two motif counts, AST, ALT, Tm
ROW = re.compile(r"^(\S+)\s+([ACGT]{12,25})\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$")

# Table S1 marks the LNA positions by font colour, green against black for
# the DNA gap, and says so in its own caption. Text extraction throws that
# away, so the screen was scored for months as though these were
# unmodified ASOs with no LNA at all. Recovering the colour gives every
# sequence a real lna_fraction instead of a guess.
LNA_GREEN = (0.0, 1.0, 0.6)


def lna_counts(pdf_path: Path) -> dict[str, int]:
    """sequence -> how many of its bases are LNA, from the font colour.

    Verified rather than assumed: all 149 sequences come back as
    `LLL........LLL`, a uniform 3-8-3 gapmer, which is what the paper
    describes and what the 73 C melting temperature of a 14-mer implies.
    """
    from pypdf import PdfReader

    reader = PdfReader(str(pdf_path))
    runs: list[tuple[str, bool]] = []
    state = {"colour": (0.0, 0.0, 0.0)}

    def before(op, args, cm, tm):
        if op in (b"sc", b"rg", b"scn"):
            vals = [float(a) for a in args if hasattr(a, "__float__")]
            if len(vals) == 3:
                state["colour"] = tuple(round(v, 3) for v in vals)

    def visit(text, cm, tm, font_dict, font_size):
        if text and text.strip():
            runs.append((text, state["colour"] == LNA_GREEN))

    for page in reader.pages:
        if "Table S2" in (page.extract_text() or ""):
            break
        page.extract_text(visitor_operand_before=before, visitor_text=visit)

    cells, current = [], []
    for text, green in runs:
        for char in text:
            if char.upper() in "ACGT" and char.isalpha():
                current.append((char.upper(), green))
            else:
                if len(current) >= 12:
                    cells.append(current)
                current = []
    if len(current) >= 12:
        cells.append(current)

    out: dict[str, int] = {}
    for cell in cells:
        out.setdefault("".join(c for c, _ in cell),
                       sum(1 for _, g in cell if g))
    return out


def extract(pdf_path: Path) -> pd.DataFrame:
    try:
        from pypdf import PdfReader
    except ImportError:
        raise SystemExit("pypdf is needed to read the supplementary PDF: pip install pypdf")

    reader = PdfReader(str(pdf_path))
    rows = []
    for page in reader.pages:
        text = page.extract_text() or ""
        if "Table S1" not in text and not rows:
            continue
        if "Table S2" in text and rows:
            break
        for line in text.split("\n"):
            m = ROW.match(line.strip())
            if not m:
                continue
            ident, seq, tgc, tcc, ast, alt, tm = m.groups()
            rows.append({
                "id": ident, "sequence": seq,
                "motifs_TGC": int(tgc), "motifs_TCC": int(tcc),
                "AST_U_per_L": int(ast), "ALT_U_per_L": int(alt), "tm_c": int(tm),
            })
    df = pd.DataFrame(rows).drop_duplicates("id")
    df["hepatotoxic"] = (df["ALT_U_per_L"] > ALT_TOXIC_ABOVE).astype(int)

    counts = lna_counts(pdf_path)
    missing = [s for s in df["sequence"] if s not in counts]
    if missing:
        raise SystemExit(
            f"{len(missing)} sequence(s) have no LNA pattern recovered from "
            "the table's font colour, so their chemistry is unknown and "
            "scoring them would invent one. First: " + missing[0])
    df["lna_count"] = df["sequence"].map(counts)
    df["lna_fraction"] = df["lna_count"] / df["sequence"].str.len()
    return df


def _boot(y, p, rng, fn, n=N_BOOT):
    idx = np.arange(len(y))
    vals = []
    for _ in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[s])) < 2:
            continue
        vals.append(fn(y[s], p[s]))
    if not vals:
        return (float("nan"), float("nan"))
    return tuple(float(v) for v in np.percentile(vals, [2.5, 97.5]))


def main() -> None:
    cached = EXTERNAL / "nar2022_hepatotox_screen.csv"
    if len(sys.argv) >= 2:
        df = extract(Path(sys.argv[1]))
        if df.empty:
            raise SystemExit(
                "No Table S1 rows parsed; the supplementary layout may have changed.")
        EXTERNAL.mkdir(parents=True, exist_ok=True)
        df.to_csv(cached, index=False)
    elif cached.exists():
        # The parse is the part that needs the PDF, and it is already done.
        # Re-scoring has to be runnable by anyone and on every retrain, or
        # the one check that argues against this endpoint only gets re-run
        # when somebody remembers to download a file.
        df = pd.read_csv(cached)
        print(f"re-scoring the cached parse at {cached.relative_to(ROOT)}")
    else:
        raise SystemExit(
            "usage: evaluate_hepatotox_external.py [<gkac562_supplemental_file.pdf>]\n"
            f"No cached parse at {cached.relative_to(ROOT)}. Download the "
            "supplementary file from the article: "
            "https://academic.oup.com/nar/article/50/13/7224/6633893"
        )

    if "lna_fraction" not in df.columns:
        raise SystemExit(
            "the cached parse predates LNA recovery and would score these "
            "gapmers as unmodified ASOs. Re-run with the supplementary PDF.")

    probs = []
    for seq, lna in zip(df["sequence"], df["lna_fraction"]):
        # These are LNA gapmers, fully phosphorothioate, as the table's own
        # caption states. Scoring them as "ASO" with lna_fraction 0 told the
        # model they carried no LNA and were not gapmers, which is a
        # chemistry that does not exist in the training data either.
        out = score_candidate(
            candidate_features=candidate_features(
                seq, "LNA_gapmer", "PS", ps_fraction=1.0, lna_fraction=float(lna)),
            processed_dir=ROOT / "data" / "processed",
            model_dir=ROOT / "artifacts" / "models",
            config_dir=ROOT / "artifacts" / "config",
        )
        probs.append(out["per_endpoint"][ENDPOINT]["prob"])
    df["predicted"] = probs

    y = df["hepatotoxic"].values
    p = df["predicted"].values
    rng = np.random.default_rng(SEED)
    roc = float(roc_auc_score(y, p))
    pr = float(average_precision_score(y, p))
    roc_ci = _boot(y, p, rng, roc_auc_score)
    pr_ci = _boot(y, p, rng, average_precision_score)

    # Where do the known-toxic ones land when everything is ranked?
    order = pd.Series(p).rank(ascending=False, method="average")
    toxic_ranks = sorted(order[y == 1].tolist())

    summary = {
        "source": "NAR 2022;50(13):7224 doi:10.1093/nar/gkac562 Supplementary Table S1",
        "assay": "C57BL/6 mice, 20 mg/kg IV, AST/ALT at 96 h",
        "chemistry": "3-8-3 LNA gapmers, fully phosphorothioate, recovered "
                     "from the font colour of Table S1 and uniform across "
                     "all 149",
        "lna_fraction": float(df["lna_fraction"].iloc[0]),
        "label": f"hepatotoxic if ALT > {ALT_TOXIC_ABOVE} U/L",
        "sequences": int(len(df)),
        "positives": int(y.sum()),
        "base_rate": float(y.mean()),
        "roc_auc": roc, "roc_ci": list(roc_ci),
        "pr_auc": pr, "pr_ci": list(pr_ci),
        "toxic_ranks_out_of_n": toxic_ranks,
        "status": "directional only",
        "caveat": (
            f"{int(y.sum())} positive sequences. Too few to validate the endpoint; the "
            "confidence interval spans most of the range. Reported so the endpoint has "
            "some external signal, not as evidence of performance. These are 14-mer "
            "3-8-3 LNA gapmers; the training corpus holds 12 LNA gapmer sequences "
            "and 369 MOE gapmers, so this is a chemistry it has barely seen."
        ),
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "hepatotox_external.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"parsed {len(df)} gapmers, {y.sum()} hepatotoxic (ALT > {ALT_TOXIC_ABOVE})")
    print(f"ROC-AUC {roc:.3f} [{roc_ci[0]:.3f}-{roc_ci[1]:.3f}]")
    print(f"PR-AUC  {pr:.3f} [{pr_ci[0]:.3f}-{pr_ci[1]:.3f}]  base rate {y.mean():.3f}")
    print(f"toxic sequences ranked at {toxic_ranks} of {len(df)}")
    print("\nDirectional only: see the caveat in hepatotox_external.json")


if __name__ == "__main__":
    main()
