import re
from pathlib import Path

print("="*80)
print("FIXING OLIGOTOX DOCUMENTATION")
print("="*80)

# The actual endpoint list from training
ACTUAL_ENDPOINTS = [
    'cns_toxicity',
    'coagulation_effect',
    'complement_activation',
    'cytokine_release',
    'hepatic_vacuolation',
    'hepatocellular_injury',
    'injection_site_reaction',
    'innate_immune_activation',
    'kupffer_cell_activation',
    'offtarget_gene_silencing',
    'renal_accumulation',
    'renal_tubular_injury',
    'thrombocytopenia'
]

ACTUAL_COUNT = len(ACTUAL_ENDPOINTS)
print(f"\nActual endpoint count: {ACTUAL_COUNT}")
print(f"Endpoints: {', '.join(ACTUAL_ENDPOINTS[:5])}... (showing first 5)")

# Files to update
files_to_update = [
    'docs/TRAINING_PIPELINE_ARCHITECTURE.md',
    'README.md'
]

for filepath in files_to_update:
    path = Path(filepath)
    if not path.exists():
        print(f"\n⚠️  {filepath} not found, skipping")
        continue
    
    print(f"\n📝 Processing {filepath}...")
    content = path.read_text()
    original_content = content
    
    # Replace "19 endpoints" with "13 endpoints"
    content = re.sub(r'19\s+(?:independent\s+)?(?:binary\s+)?endpoints?', 
                     f'{ACTUAL_COUNT} endpoints', content, flags=re.IGNORECASE)
    content = re.sub(r'19\s+toxicity\s+endpoints?', 
                     f'{ACTUAL_COUNT} toxicity endpoints', content, flags=re.IGNORECASE)
    
    # Replace "19 independent binary classifiers"
    content = re.sub(r'\*\*19\s+independent\s+binary\s+classifiers\*\*',
                     f'**{ACTUAL_COUNT} independent binary classifiers**', content)
    
    if content != original_content:
        # Create backup
        backup_path = path.with_suffix('.md.bak')
        backup_path.write_text(original_content)
        print(f"  ✓ Created backup: {backup_path}")
        
        # Write updated content
        path.write_text(content)
        print(f"  ✓ Updated {filepath}")
        print(f"    Changed '19 endpoints' → '{ACTUAL_COUNT} endpoints'")
    else:
        print(f"  ℹ️  No changes needed in {filepath}")

print("\n" + "="*80)
print("DOCUMENTATION UPDATE COMPLETE")
print("="*80)
print(f"\n✅ All files processed")
print(f"📊 Platform now correctly documents {ACTUAL_COUNT} endpoints")
