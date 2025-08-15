# Fully Decarbonising the Off-Grid Al Jalamid Phosphate Complex
**MILP Formulation and Techno-Economics for PV–CSP–BESS–PEM Integration**

This repository contains the models, data, and documentation supporting the MSc Business Analytics Individual Research Report (IRR) on optimising a hybrid renewable–hydrogen system for Ma’aden's off-grid Al Jalamid phosphate complex.

## 📜 Overview
We develop a **two-pass lexicographic Mixed-Integer Linear Program (MILP)** to size and dispatch:
- **Photovoltaic (PV)** field
- **Concentrating Solar Power (CSP)** with molten-salt Thermal Energy Storage (TES)
- **Battery Energy Storage System (BESS)**
- **Proton-Exchange-Membrane (PEM) Electrolyser**

The system replaces diesel-based electricity and Heavy Fuel Oil (HFO) process heat with renewable electricity and green hydrogen.

## 📊 Key Features
- **Deterministic and CVaR (risk-aware) MILP** formulations
- Hourly resolution over 1-year and 10-year NASA POWER irradiance records
- Explicit coupling of electrolyser hydrogen output to avoided HFO use and carbon credits
- Sensitivity analysis for **CO₂ price** and **production scheduling**
- Open, auditable Excel outputs (hourly, monthly, summary)

## 📂 Repository Structure
```
📁 data/           → Input datasets (irradiance, load, techno-economic parameters)
📁 src/            → MATLAB MILP formulation and drivers
📁 results/        → Excel outputs for deterministic and CVaR runs
📁 docs/           → IRR PDF and supplementary appendices
README.md          → This file
LICENSE            → License for reuse
```

## 🖥 Requirements
- MATLAB R2023a+ with Optimization Toolbox
- NASA POWER API access (optional, for data refresh)

## 🚀 Running the Model
1. Place `.mat` irradiance and operational datasets in `data/`
2. Adjust parameters in `unit` and `econ` structs
3. Run:
```matlab
run_milp_milp_v2_lex.m
```
4. Outputs are stored in `results/` in Excel format.

## 📄 License
MIT License – see [LICENSE](LICENSE) file.

## 📚 Citation
If you use this code or data, please cite:
> Moreira, F. (2025). *Fully Decarbonising the Off-Grid Al Jalamid Phosphate Complex: MILP Formulation and Techno-Economics*. MSc IRR, Imperial College London.
