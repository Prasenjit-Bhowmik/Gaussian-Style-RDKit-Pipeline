# -----------------------------
# mini_gaussian_pipeline.py
# Stage 1 → 6 (Full Gaussian-style pipeline, RDKit-safe)
# optional Stage 7 QM
# -----------------------------

import os
import pandas as pd
import shutil

def safe_name(name):
    return name.replace(" ", "_").replace("(", "").replace(")", "")

def xtb_available():
    return shutil.which("xtb") is not None

def orca_available():
    return shutil.which("orca") is not None

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from mordred import Calculator, descriptors
from pathlib import Path

# -----------------------------
# PARAMETERS
# -----------------------------
smiles_dir = Path("smiles")        # folder with .smi files
mol2d_dir = smiles_dir / "mol_files_2D"
mol3d_dir = smiles_dir / "mol_files_3D"
n_conformers = 20                  # Stage 3: number of conformers

desc2d_csv = smiles_dir / "all_2D_descriptors.csv"
desc3d_csv = smiles_dir / "all_3D_descriptors.csv"
prop_csv = smiles_dir / "all_properties.csv"

# Create folders
mol2d_dir.mkdir(exist_ok=True)
mol3d_dir.mkdir(exist_ok=True)

# -----------------------------
# Stage 1 + 2 + 3: SMILES → 2D + 3D + lowest-energy
# -----------------------------
smi_files = [f for f in os.listdir(smiles_dir) if f.endswith(".smi")]
mols_2d = []
names = []

mols_3d = []
energies = []

for smi_file in smi_files:
    with open(smiles_dir / smi_file) as f:
        smi = f.read().strip()

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"⚠ Failed: {smi_file}")
        continue

    name = smi_file.replace(".smi", "")
    names.append(name)
    mols_2d.append(mol)

    # 2D MOL
    Chem.MolToMolFile(mol, str(mol2d_dir / f"{name}.mol"))

    # 3D + multiple conformers
    mol3d = Chem.AddHs(mol)
    conf_ids = AllChem.EmbedMultipleConfs(mol3d, numConfs=n_conformers, params=AllChem.ETKDGv3())

    conf_energies = []
    for conf_id in conf_ids:
        if AllChem.MMFFHasAllMoleculeParams(mol3d):
            ff = AllChem.MMFFGetMoleculeForceField(mol3d, AllChem.MMFFGetMoleculeProperties(mol3d), confId=conf_id)
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol3d, confId=conf_id)
        energy = ff.CalcEnergy()
        ff.Minimize()
        conf_energies.append((conf_id, energy))

    # Lowest-energy conformer
    lowest_conf_id, lowest_energy = min(conf_energies, key=lambda x: x[1])
    energies.append(lowest_energy)
    Chem.MolToMolFile(mol3d, str(mol3d_dir / f"{name}_lowest.mol"), confId=lowest_conf_id)
    mols_3d.append(mol3d)

print(f"✅ {len(mols_2d)} molecules processed: 2D + 3D + lowest-energy conformers")

# -----------------------------
# Stage 4: 2D descriptors
# -----------------------------
calc2d = Calculator(descriptors, ignore_3D=True)
df2d = calc2d.pandas(mols_2d)
df2d.insert(0, "Molecule", names)
df2d.to_csv(str(desc2d_csv), index=False)
print(f"✅ 2D descriptors saved: {desc2d_csv}")

# -----------------------------
# Stage 4: 3D descriptors + energies
# -----------------------------
calc3d = Calculator(descriptors, ignore_3D=False)
df3d = calc3d.pandas(mols_3d)
df3d.insert(0, "Molecule", names)
df3d["Energy_kcal_mol"] = energies
df3d.to_csv(str(desc3d_csv), index=False)
print(f"✅ 3D descriptors + energies saved: {desc3d_csv}")

# -----------------------------
# Stage 5: Physicochemical properties
# -----------------------------
properties = []
for mol, name, energy in zip(mols_3d, names, energies):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    lipinski = (mw <= 500) and (logp <=5) and (hbd <=5) and (hba <=10)

    properties.append({
        "Molecule": name,
        "MW": mw,
        "LogP": logp,
        "TPSA": tpsa,
        "HBD": hbd,
        "HBA": hba,
        "Lipinski": lipinski,
        "Energy_kcal_mol": energy
    })

df_prop = pd.DataFrame(properties)
df_prop.to_csv(str(prop_csv), index=False)
print(f"✅ Physicochemical properties saved: {prop_csv}")

# -----------------------------
# Stage 6: Batch + one-command execution
# -----------------------------
print("✅ Pipeline completed: 2D MOLs, 3D MOLs, descriptors, energies, properties")
print("All outputs saved in the 'smiles/' folder. Run on any folder of .smi files with one command.")

# -----------------------------
# Stage 7: Cloud / HPC QM (xTB or ORCA) optional
# -----------------------------
print("\n🔬 Cloud/HPC QM stage starting...")

if not xtb_available():
    print("⚠ xTB not available in this environment. Skipping xTB optimization.")
else:
    print("✅ xTB found. Running GFN2-xTB optimizations for all molecules...")
    for name in names:
        sname = safe_name(name)
        xtb_dir = smiles_dir / "qm_xtb"
        work = xtb_dir / sname
        work.mkdir(parents=True, exist_ok=True)

        molfile = mol3d_dir / f"{name}_lowest.mol"
        cmd = f'xtb "{molfile}" --gfn 2 --opt > xtb.out'
        os.system(f'cd "{work}" && {cmd}')

    print("✅ xTB optimization completed. Output stored in smiles/qm_xtb/")

# Optional ORCA stage
if orca_available():
    run_orca = input("Run ORCA DFT on shortlisted molecules? (yes/no): ").strip().lower()
    if run_orca == "yes":
        print("✅ ORCA found. Preparing ORCA jobs...")
        # Here you can add your ORCA input generation and submission logic
        print("⚠ ORCA execution logic not implemented yet. Add HPC submission code here.")
else:
    print("⚠ ORCA not available in this environment. Skipping ORCA stage.")

print("🔬 Cloud/HPC QM stage finished.\n")





