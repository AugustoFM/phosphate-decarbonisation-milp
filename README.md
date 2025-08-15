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
- 📁 data/ → Input datasets (irradiance, load, techno-economic parameters)
- 📁 src/ → MATLAB MILP formulation and drivers
- 📁 results/ → Excel outputs for deterministic and CVaR runs
- 📁 docs/ → IRR PDF and supplementary appendices
