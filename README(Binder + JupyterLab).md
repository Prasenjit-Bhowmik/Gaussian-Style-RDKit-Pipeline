
---

## ▶ Run in Binder (One-Click)

[Open and run pipeline in Binder](https://mybinder.org/v2/gh/Prasenjit-Bhowmik/mini_gaussian_pipeline/main?filepath=run_pipeline.ipynb)

1. Binder will launch JupyterLab in the cloud.  
2. Open `run_pipeline.ipynb` notebook.  
3. Click **Run All** → pipeline runs automatically.  
4. All outputs are saved in the `smiles/` folder.  

---

## 💻 Requirements

- Python ≥ 3.11  
- Packages: `rdkit-pypi`, `mordred`, `pandas`, `numpy`  
- All dependencies are installed automatically by Binder.  

---

## ⚡ Quick Local Run

If you want to run locally:

```bash
git clone https://github.com/Prasenjit-Bhowmik/mini_gaussian_pipeline.git
cd mini_gaussian_pipeline
pip install -r requirements.txt
python mini_gaussian_pipeline.py
