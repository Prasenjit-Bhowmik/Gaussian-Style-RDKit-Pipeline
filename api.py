from fastapi import FastAPI
from xtb_runner import run_xtb

app = FastAPI()

@app.post("/optimize")
def optimize(data: dict):
    xyz, energy = run_xtb(
        molfile=data["molfile"],
        charge=data.get("charge", 0)
    )
    return {
        "energy_Eh": energy,
        "optimized_xyz": xyz
    }
