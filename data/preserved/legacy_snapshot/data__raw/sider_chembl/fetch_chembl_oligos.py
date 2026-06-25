"""Fetch ChEMBL records for approved oligonucleotide drugs."""
from chembl_webresource_client.new_client import new_client
import json, os

DRUG_NAMES = [
    "patisiran","vutrisiran","givosiran","lumasiran","inclisiran",
    "nedosiran","fitusiran","plozasiran","inotersen","eplontersen",
    "tofersen","nusinersen","mipomersen","olezarsen","donidalorsen",
    "eteplirsen","golodirsen","viltolarsen","casimersen","pegaptanib","avacincaptad"
]
molecule = new_client.molecule
out = []
for name in DRUG_NAMES:
    results = molecule.filter(pref_name__iexact=name).only(
        ['molecule_chembl_id','pref_name','molecule_type',
         'max_phase','black_box_warning','oral','parenteral',
         'molecule_properties'])
    r = list(results)
    if r:
        out.append(r[0])
        print(f"  ✓ {name}: {r[0].get('molecule_chembl_id')}")
    else:
        print(f"  - {name}: not found in ChEMBL")

out_path = os.path.join(os.path.dirname(__file__), "chembl_approved_oligos.json")
with open(out_path, "w") as f:
    json.dump(out, f, indent=2)
print(f"\nSaved {len(out)} records to {out_path}")
