import gradio as gr
import os
from pathlib import Path
from mini_gaussian_pipeline import run_pipeline
import shutil
import zipfile

# Temporary folder for user uploads
USER_FOLDER = Path("smiles")
USER_FOLDER.mkdir(parents=True, exist_ok=True)

def run_gradio_pipeline(smi_file, smi_text):
    # Clear previous inputs
    if USER_FOLDER.exists():
        shutil.rmtree(USER_FOLDER)
    USER_FOLDER.mkdir(parents=True, exist_ok=True)

    # Save uploaded file if provided
    if smi_file is not None:
        file_path = USER_FOLDER / smi_file.name
        smi_file.save(file_path)

    # Save pasted text as input.smi if provided
    if smi_text:
        file_path = USER_FOLDER / "input.smi"
        with open(file_path, "w") as f:
            f.write(smi_text)

    # Run your pipeline
    run_pipeline(str(USER_FOLDER))

    # Create a zip of all output CSVs
    zip_path = USER_FOLDER / "results.zip"
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        for csv_file in ["all_2D_descriptors.csv", "all_3D_descriptors.csv", "all_properties.csv"]:
            fpath = USER_FOLDER / csv_file
            if fpath.exists():
                zipf.write(fpath, arcname=csv_file)

    return zip_path

# Gradio interface
iface = gr.Interface(
    fn=run_gradio_pipeline,
    inputs=[
        gr.File(label="Upload .smi file", file_types=['.smi']),
        gr.Textbox(label="Or paste SMILES here (one per line)")
    ],
    outputs=gr.File(label="Download results (ZIP)"),
    title="Gaussian-style RDKit Pipeline",
    description="Upload a .smi file or paste SMILES. The pipeline generates 2D/3D structures, descriptors, and physicochemical properties.",
    allow_flagging="never"
)

if __name__ == "__main__":
    iface.launch()
