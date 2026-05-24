# Academic Epidemiology Dashboard v5.1
## Zoonotic Disease Surveillance Portal (ZDSP)

---

> **Authorship Notice**
>
> This software was independently conceived and developed by **Felix Loaiza** as a personal hobby project, not commissioned or published by any academic institution, government agency, or public health organization. It is shared openly as a portfolio demonstration of applied competency in computational epidemiology and virology.
>
> All epidemiological outputs, projections, and simulation results are **synthetically generated or AI-inferred from historical, publicly available datasets**. No real-time surveillance feeds or live clinical data are connected to this system.

---

## Overview

A research-inspired web application for the exploratory analysis and visualization of zoonotic disease transmission dynamics. Built with Python (Dash / Plotly), it integrates stochastic SEIR modeling, genomic surveillance interpretation, wastewater-based epidemiology (WBE) simulation, and environmental forcing analysis into a unified dark academic interface.

The diseases implemented cover distinct transmission archetypes: rodent-borne spillover (ANDV) and direct-contact hemorrhagic fever (EBOV). This range allows the modeling framework to be tested across epidemiologically diverse scenarios.

---

## Features

### Stochastic SEIR Projections
180-day forward projections via Monte Carlo sampling (500 to 50,000 samples). Renders 95% and 50% credible interval bands following CDC/WHO visualization conventions.

### Endemic Channel Alert System (Bortman Method)
Seasonal baselines across 52 epidemiological weeks. Four risk zones derived from historical percentile thresholds:

| Zone | Threshold | Meaning |
|---|---|---|
| Success | < P25 | Below-average incidence |
| Security | P25 to P50 | Expected baseline |
| Alert | P50 to P75 | Early warning |
| Epidemic | > P75 | Emergency threshold |

### Genomic Surveillance Suite
Illustrative NCBI GenBank accession IDs with phylogenetic clade assignment. dN/dS ratio monitoring for positive selection detection. Nucleotide diversity (pi) tracking across synthetic variant datasets.

### Wastewater-Based Epidemiology (WBE)
RT-qPCR and RT-ddPCR assay simulation with PMMoV-normalized genomic copy estimates. Simulates early signal detection relative to hypothetical clinical reporting thresholds.

### Environmental Forcing Index (EFI)
Temperature, humidity, and precipitation correlation analysis sourced via Open-Meteo API (with ERA5 reanalysis fallback). Visualized as EFI vs. effective reproduction number (Rt) scatter plot.

### Disease Parameter Suite
Spillover rate (beta_z), human-to-human transmission rate (beta_h), recovery rate (gamma), CFR from peer-reviewed literature, and R0 with risk-stratified color coding.

---

## Supported Diseases

### Hantavirus Andes (ANDV)

| Parameter | Value |
|---|---|
| Classification | Zoonotic / Rodent Contact |
| Scientific name | Andes orthohantavirus |
| Primary reservoir | Oligoryzomys longicaudatus |
| Model type | Zoonotic SEIR |
| beta_z | 0.0500 |
| beta_h | 0.0050 |
| gamma | 0.1667 |
| CFR | 35% (Ferres et al., 2007) |
| Geographic focus | Chile, Argentina |
| Peak season | Austral summer, weeks 1 to 10 and 45 to 52 |
| dN/dS | 0.43 (Purifying selection) |
| GenBank reference | NC_003466.1 (illustrative) |
| Climate data | ERA5 normals, Valdivia (-41.14, -72.78) |

### Ebola Virus Disease (BDBV / EBOV)

| Parameter | Value |
|---|---|
| Classification | Zoonotic / Direct Contact |
| Scientific name | Zaire ebolavirus |
| Model type | Zoonotic SEIRD (with death compartment) |
| beta_z | 0.0100 |
| beta_h | 0.1250 |
| gamma | 0.0714 |
| CFR | 57.5% (WHO validated range) |
| R0 (human) | 1.75 (Althaus, 2014; empirical range 1.5 to 2.5) |
| Geographic focus | DRC, West Africa |
| GenBank reference | KM034562.1, EBOV-Makona lineage (illustrative) |
| dN/dS | 0.18 (Strong purifying selection) |
| EFI score | 0.89 (Very High Risk) |
| Climate data | ERA5 normals, Kinshasa (-4.32, 15.32) |

---

## Architecture

### Files

| File | Role |
|---|---|
| `academic_dashboard_app.py` | Dash application, layout, callbacks, visualization functions |
| `scientific_epidemiology_backend.py` | SEIR solver, genomic module, climate module, WBE module, endemic channel calculator |
| `generate_historical_csvs.py` | Generates synthetic historical incidence CSVs for the Bortman endemic channel |
| `data_cache/` | Output directory for generated CSVs and cached API responses |

### Backend Modules

**GenomicSurveillanceModule** -- Fetches variant metadata via NCBI Entrez API. Calculates nucleotide diversity (pi) and dN/dS ratios. Falls back to registry values when offline.

**EcologicalSurveillanceModule** -- Retrieves climate data from Open-Meteo API with exponential backoff (tenacity). Falls back to ERA5 reanalysis monthly normals per disease geographic focus. Computes the Environmental Forcing Index.

**WastewaterSurveillanceModule** -- Simulates RT-qPCR / RT-ddPCR signal with PMMoV normalization. Generates WBE reports per disease key.

**ZoonoticSEIRModel** -- Stochastic ODE solver for the zoonotic SEIR/SEIRD system. Supports Monte Carlo sampling for credible interval generation.

**calculate_endemic_channels()** -- Implements the Bortman method over historical weekly incidence CSVs. Returns P25/P50/P75 percentile bands per epidemiological week.

### Dashboard Visualization Functions

| Function | Output |
|---|---|
| `create_confidence_interval_chart()` | 95% and 50% CrI stochastic projection bands |
| `create_endemic_channel_chart()` | Bortman seasonal baseline with four risk zones |
| `create_efi_rt_chart()` | EFI vs. Rt scatter plot |
| `create_dnds_alert_panel()` | Genomic selection pressure monitor |
| `create_genomic_surveillance_table()` | Variant registry table |

### Mathematical Model

```
dS/dt = -beta_z * Z(t) * S/N - beta_h * I * S/N
dE/dt = +beta_z * Z(t) * S/N + beta_h * I * S/N - sigma * E
dI/dt = +sigma * E - gamma * I
dR/dt = +gamma * I
```

R0 = beta_h / gamma

| R0 | Risk |
|---|---|
| < 1.0 | Elimination trajectory |
| 1.0 to 2.5 | Endemic baseline |
| 2.5 to 5.0 | Serious outbreak risk |
| > 5.0 | Critical pandemic threshold |

dN/dS interpretation: values well below 1 indicate purifying selection (functional constraint); values approaching or exceeding 1 signal adaptive evolution and elevated pandemic risk.

---

## Historical CSV Generation

`generate_historical_csvs.py` produces synthetic weekly incidence series (cases/100k/week) for the Bortman endemic channel. Only Hantavirus Andes and Ebola EVD have full CSV pipelines. Dengue DENV-2 and COVID-19 use hardcoded synthetic seasonal curves defined directly in the dashboard.

**Hantavirus Andes:** Based on MINSAL Chile bulletins 2015 to 2024. Annual incidence of 300 to 600 cases in Chile (~19M population). Lognormal noise plus stochastic outbreak events (2.5% of weeks, 3x to 5x baseline). Seed: 20260523.

**Ebola EVD:** Based on WHO Situation Reports 2018 to 2020 (North Kivu / Ituri outbreak). Documented outbreak windows simulated for 2018, 2019, 2020, 2021, and 2022 using a sine-shaped epidemic pulse. Interdemic baseline ~0.0015 cases/100k/week.

---

## Dependencies

```bash
pip install dash plotly numpy pandas requests tenacity pymc
```

| Library | Role |
|---|---|
| dash >= 2.0 | Web framework and reactive callbacks |
| plotly >= 5.0 | Interactive charts |
| numpy >= 1.20 | Numerical computation and Monte Carlo sampling |
| pandas >= 1.3 | Data manipulation |
| requests | NCBI Entrez and Open-Meteo API calls |
| tenacity | Exponential backoff retry policy for external APIs |
| pymc | Bayesian parameter estimation (optional) |

---

## Running the Application

```bash
python generate_historical_csvs.py   # Generate historical incidence CSVs
python scientific_epidemiology_backend.py     # Generate historical surveillance report
python academic_dashboard_app.py     # Launch dashboard
```

Available at `http://localhost:8050`

The dashboard falls back to deterministic synthetic curves if `scientific_epidemiology_backend.py` is absent. With it present, live SEIR simulations, climate API calls, and genomic module outputs are enabled.

---

## Color Coding

| Status | Hex | Meaning |
|---|---|---|
| Normal | #10b981 | Purifying selection / Low R0 / Below baseline |
| Alert | #f59e0b | Watch zone / Elevated R0 / Early warning |
| Critical | #ef4444 | Pandemic threshold / High CFR / Epidemic phase |
| Info | #3b82f6 | Neutral / Baseline / Security zone |

---

## Data Provenance and Limitations

All simulation outputs are synthetic. Sources:

- AI-assisted inference from historical, publicly available epidemiological literature
- Stochastic and deterministic mathematical models parameterized from peer-reviewed publications
- Open-Meteo API climate data and ERA5 reanalysis fallback normals
- NCBI Entrez API for illustrative genomic accession references

No real-time surveillance systems or clinical databases are connected.

| Limitation | Detail |
|---|---|
| Synthetic projections | Deterministic fallback curves when backend is absent |
| MCMC sample size | Default 500 samples; production standard is 5,000 to 50,000 |
| No real-time nowcasting | Not calibrated for operational outbreak response |
| CFR point estimates | Literature-derived; real values vary by population and healthcare access |
| WBE module | Simulation only; requires lab calibration for operational use |
| AI-inferred parameters | Some values derived from LLM inference on historical datasets |

---

## References

- Bortman, M. (1999). Elaboracion de corredores o canales endemicos mediante planillas de calculo. Revista Panamericana de Salud Publica, 5(1), 1-8.
- Althaus, C. L. (2014). Estimating the reproduction number of Ebola virus (EBOV) during the 2014 outbreak in West Africa. PLOS Currents Outbreaks.
- WHO Ebola Response Team. (2014). Ebola virus disease in West Africa -- the first 9 months. New England Journal of Medicine, 371, 1481-1495.
- Ferres, M., et al. (2007). Prospective evaluation of household contacts of persons with hantavirus cardiopulmonary syndrome in Chile. New England Journal of Medicine, 357, 2564-2571.
- Hersbach, H., et al. (2020). The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146, 1999-2049.
- NCBI GenBank: https://www.ncbi.nlm.nih.gov/genbank/

---

## Changelog

### v5.1
- Full English translation of codebase
- Ebola Virus Disease integration (BDBV/EBOV)
- 95% and 50% credible interval bands
- Bortman endemic channel with historical CSVs
- dN/dS genomic selection monitoring
- Environmental Forcing Index with Open-Meteo integration and ERA5 fallback
- WBE simulation module

### v5.0
- Backend architecture refactoring
- SEIR parameter optimization

### v4.x
- Initial release with Dengue and Hantavirus support

---

## Disclaimer

Independent hobby project by **Felix Loaiza**. Not affiliated with any academic institution, research organization, or public health authority. All outputs are synthetic and must not be used for clinical decision-making or operational epidemiological surveillance.

**Author:** Felix Loaiza | **Version:** 1.0 | **Released:** May 2026
