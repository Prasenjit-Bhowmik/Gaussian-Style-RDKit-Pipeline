import gradio as gr
from pathlib import Path
import shutil
import pandas as pd
from mini_gaussian_pipeline import mini_gaussian_pipeline  # you wrap your main function

# Temporary working folder
WORK_DIR = Path("temp_work")
WORK_DIR.mkdir(exist_ok=True)

def run_pipeline(smi_file, smi_text):
    # Clear previous run
    if WORK_DIR.exists():
        shutil.rmtree(WORK_DIR)
    WORK_DIR.mkdir(exist_ok=True)

    # Save SMILES input
    if smi_file:
        content = smi_file.read().decode("utf-8")
    elif smi_text:
        content = smi_text
    else:
        return "⚠ No input provided."

    smi_path = WORK_DIR / "input.smi"
    smi_path.write_text(content)

    # Run your existing pipeline (Stages 1–6)
    mini_gaussian_pipeline(str(smi_path), output_dir=str(WORK_DIR))

    # Prepare results (ZIP)
    import zipfile
    zip_path = WORK_DIR / "results.zip"
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        for f in WORK_DIR.glob("*.csv"):
            zipf.write(f, f.name)
        for f in WORK_DIR.glob("mol_files_2D/*.mol"):
            zipf.write(f, f"mol_files_2D/{f.name}")
        for f in WORK_DIR.glob("mol_files_3D/*.mol"):
            zipf.write(f, f"mol_files_3D/{f.name}")

    return zip_path

# Gradio UI
iface = gr.Interface(
    fn=run_pipeline,
    inputs=[
        gr.File(label="Upload .smi file", file_types=[".smi"]),
        gr.Textbox(lines=5, label="Or paste SMILES (one per line)")
    ],
    outputs=gr.File(label="Download results (ZIP)"),
    title="Gaussian-style RDKit Pipeline",
    description="Upload a .smi file or paste SMILES. Generates 2D/3D MOLs, descriptors, energies, and properties."
)

if __name__ == "__main__":
    iface.launch()
