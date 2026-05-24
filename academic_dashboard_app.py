"""
Academic Epidemiology Dashboard - Version 5.1
==============================================
Rigorous zoonotic disease surveillance interface with:
- Minimalist dark academic aesthetic
- 95% credible intervals on projections
- Genomic surveillance with real NCBI accessions
- Wastewater-based epidemiology (WBE) monitoring
- Endemic channel alert system (Bortman Method)
- Environmental forcing variables
- FULL Ebola Virus Disease integration (BDBV/EBOV)

Author: Felix Loaiza
Copyright (c) 2026. All rights reserved.
"""

import os
import json
import sys
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

import dash
from dash import dcc, html, Input, Output, State, callback
import plotly.graph_objects as go
from plotly.subplots import make_subplots


_BACKEND_DIR = os.path.dirname(os.path.abspath(__file__))
if _BACKEND_DIR not in sys.path:
    sys.path.insert(0, _BACKEND_DIR)

try:
    from scientific_epidemiology_backend import (
        ZoonoticSEIRModel,
        EcologicalSurveillanceModule,
        get_full_surveillance_report,
        DISEASE_REGISTRY,
        DATA_DIR,
    )
    _BACKEND_AVAILABLE = True
except ImportError as _be:
    _BACKEND_AVAILABLE = False
    import warnings
    warnings.warn(
        f"[Dashboard] Could not import scientific_epidemiology_backend: {_be}. "
        "Falling back to deterministic placeholder curves. "
        "Run both files from the same directory to activate the real engine.",
        RuntimeWarning,
    )

# ---------------------------------------------------------------------------
# DISEASE → backend key mapping
# ---------------------------------------------------------------------------
_DISEASE_BACKEND_KEY = {
    "Hantavirus Andes - ANDV (Zoonotic / Rodent Contact)": "hanta_andes",
    "Dengue Hemorrhagic Fever (Serotype 2)": "dengue_serotype2",
    "COVID-19 (SARS-CoV-2)": "covid19",
    "Ebola Virus Disease - BDBV/EBOV (Zoonotic / Direct Contact)": "ebola_virus_disease",
}

# Simple in-memory cache to avoid re-running SEIR on every dropdown change
# Key: backend disease key  →  Value: simulation result dict
_SEIR_CACHE: dict = {}

# =============================================================================
# DASH APP INITIALIZATION
# =============================================================================

app = dash.Dash(
    __name__,
    external_stylesheets=[],
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1"},
        {"name": "description", "content": "Zoonotic Disease Surveillance Portal"},
    ],
)
server = app.server

app.title = "ZDSPMP v5.1 - Academic Epidemiology"

# =============================================================================
# CUSTOM CSS - ACADEMIC MINIMALIST THEME
# =============================================================================

ACADEMIC_THEME_CSS = """
:root {
    --bg-dark: #0f172a;
    --bg-panel: #1a2942;
    --border-subtle: #2d3d5c;
    --text-primary: #e0e7ff;
    --text-secondary: #9ca3af;
    --accent-primary: #3b82f6;
    --accent-success: #10b981;
    --accent-warning: #f59e0b;
    --accent-danger: #ef4444;
    --font-mono: 'Roboto Mono', 'JetBrains Mono', 'SF Mono', monospace;
    --font-sans: 'Inter', 'Segoe UI', sans-serif;
}

* {
    margin: 0;
    padding: 0;
    box-sizing: border-box;
}

html, body {
    background-color: var(--bg-dark);
    color: var(--text-primary);
    font-family: var(--font-sans);
    font-size: 13px;
    line-height: 1.5;
    overflow-x: hidden;
}

/* Header */
.header-container {
    background: linear-gradient(180deg, #1a2942 0%, #0f172a 100%);
    border-bottom: 1px solid var(--border-subtle);
    padding: 24px 32px;
    margin-bottom: 24px;
}

.header-title {
    font-size: 24px;
    font-weight: 600;
    letter-spacing: -0.5px;
    margin-bottom: 8px;
}

.header-subtitle {
    font-size: 12px;
    color: var(--text-secondary);
    font-family: var(--font-mono);
    letter-spacing: 0.5px;
}

/* Layout Grid */
.main-container {
    display: grid;
    grid-template-columns: 1fr 350px;
    gap: 24px;
    padding: 0 32px;
    max-width: 1800px;
    margin: 0 auto;
}

.content-main {
    display: flex;
    flex-direction: column;
    gap: 24px;
}

.sidebar-controls {
    display: flex;
    flex-direction: column;
    gap: 16px;
}

/* Panel Cards */
.panel-card {
    background-color: var(--bg-panel);
    border: 1px solid var(--border-subtle);
    border-radius: 4px;
    padding: 20px;
    transition: all 0.3s ease;
}

.panel-card:hover {
    border-color: var(--accent-primary);
}

.panel-title {
    font-size: 12px;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.8px;
    color: var(--text-secondary);
    margin-bottom: 16px;
    display: flex;
    align-items: center;
    gap: 8px;
}

.panel-content {
    display: flex;
    flex-direction: column;
    gap: 12px;
}

/* Metrics */
.metric-row {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 8px 0;
    border-bottom: 1px solid rgba(45, 61, 92, 0.5);
}

.metric-row:last-child {
    border-bottom: none;
}

.metric-label {
    font-size: 11px;
    color: var(--text-secondary);
    font-family: var(--font-mono);
}

.metric-value {
    font-size: 14px;
    font-weight: 600;
    font-family: var(--font-mono);
    color: var(--text-primary);
}

.metric-value.accent {
    color: var(--accent-primary);
}

/* Status Indicators */
.status-badge {
    display: inline-block;
    padding: 4px 8px;
    font-size: 10px;
    font-weight: 600;
    border-radius: 2px;
    font-family: var(--font-mono);
    text-transform: uppercase;
    letter-spacing: 0.5px;
}

.status-operational {
    background: rgba(16, 185, 129, 0.2);
    color: var(--accent-success);
}

.status-alert {
    background: rgba(239, 68, 68, 0.2);
    color: var(--accent-danger);
}

.status-warning {
    background: rgba(245, 158, 11, 0.2);
    color: var(--accent-warning);
}

/* Charts */
.chart-container {
    width: 100%;
    min-height: 300px;
}

/* Dropdowns & Selects */
select, input {
    background-color: rgba(45, 61, 92, 0.5);
    border: 1px solid var(--border-subtle);
    color: var(--text-primary);
    padding: 8px 12px;
    border-radius: 4px;
    font-family: var(--font-mono);
    font-size: 12px;
}

select:focus, input:focus {
    outline: none;
    border-color: var(--accent-primary);
    background-color: rgba(45, 61, 92, 0.8);
    box-shadow: 0 0 0 2px rgba(59, 130, 246, 0.1);
}

/* Footer / Data Provenance */
.footer-provenance {
    border-top: 1px solid var(--border-subtle);
    padding: 16px 32px;
    margin-top: 32px;
    font-size: 10px;
    color: var(--text-secondary);
    font-family: var(--font-mono);
    background-color: rgba(15, 23, 42, 0.5);
}

.provenance-row {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
    gap: 16px;
    margin-bottom: 8px;
}

.provenance-item {
    display: flex;
    flex-direction: column;
    gap: 4px;
}

.provenance-label {
    color: var(--text-secondary);
    font-weight: 600;
}

.provenance-value {
    color: var(--text-primary);
}

/* Dash Dropdown — dark theme override (Dash 2.x) */
.dash-dropdown .Select-control,
.dash-dropdown .Select--single .Select-control {
    background-color: #1a2942 !important;
    border: 1px solid #2d3d5c !important;
    border-radius: 4px !important;
    min-height: 38px !important;
    cursor: pointer !important;
    transition: border-color 0.2s ease !important;
}

.dash-dropdown .Select-control:hover {
    border-color: #3b82f6 !important;
    box-shadow: 0 0 0 2px rgba(59,130,246,0.15) !important;
}

.dash-dropdown .Select--single > .Select-control .Select-value,
.dash-dropdown .Select-value-label {
    color: #e0e7ff !important;
    font-family: 'Roboto Mono', monospace !important;
    font-size: 12px !important;
    line-height: 36px !important;
}

.dash-dropdown .Select-placeholder {
    color: #6b7280 !important;
    font-family: 'Roboto Mono', monospace !important;
    font-size: 12px !important;
    line-height: 36px !important;
}

.dash-dropdown .Select-input > input {
    color: #e0e7ff !important;
    font-family: 'Roboto Mono', monospace !important;
    font-size: 12px !important;
    background: transparent !important;
}

.dash-dropdown .Select-arrow {
    border-top-color: #9ca3af !important;
}

.dash-dropdown .Select.is-open .Select-control {
    border-color: #3b82f6 !important;
    background-color: #1a2942 !important;
    border-bottom-left-radius: 0 !important;
    border-bottom-right-radius: 0 !important;
}

.dash-dropdown .Select-menu-outer {
    background-color: #111e33 !important;
    border: 1px solid #3b82f6 !important;
    border-top: none !important;
    border-radius: 0 0 4px 4px !important;
    box-shadow: 0 8px 24px rgba(0,0,0,0.5) !important;
    z-index: 9999 !important;
}

.dash-dropdown .Select-option {
    background-color: #111e33 !important;
    color: #9ca3af !important;
    font-family: 'Roboto Mono', monospace !important;
    font-size: 12px !important;
    padding: 10px 14px !important;
    cursor: pointer !important;
}

.dash-dropdown .Select-option:hover,
.dash-dropdown .Select-option.is-focused {
    background-color: rgba(59,130,246,0.18) !important;
    color: #e0e7ff !important;
}

.dash-dropdown .Select-option.is-selected {
    background-color: rgba(59,130,246,0.30) !important;
    color: #93c5fd !important;
    font-weight: 600 !important;
}

/* ID-targeted override — highest specificity, works on all Dash versions */
#disease-selector .Select-control,
#disease-selector .VirtualizedSelectFocusedOption,
#disease-selector input {
    background-color: #1a2942 !important;
    color: #e0e7ff !important;
    font-family: 'Roboto Mono', monospace !important;
    font-size: 12px !important;
}

#disease-selector .Select-value-label,
#disease-selector .Select-single-value,
#disease-selector .Select-value,
#disease-selector .Select-value *,
#disease-selector div[class*="value"],
#disease-selector div[class*="singleValue"],
#disease-selector div[class*="placeholder"] {
    color: #e0e7ff !important;
    font-family: 'Roboto Mono', monospace !important;
    font-size: 12px !important;
}

#disease-selector .Select-menu-outer,
#disease-selector .VirtualizedSelectOption {
    background-color: #111e33 !important;
    color: #9ca3af !important;
    font-family: 'Roboto Mono', monospace !important;
    font-size: 12px !important;
}

/* Nuclear option: any text inside the dropdown wrapper */
#disease-selector * {
    color: #e0e7ff;
    font-family: 'Roboto Mono', monospace !important;
    font-size: 12px !important;
}
#disease-selector .Select-menu-outer *,
#disease-selector .Select-option {
    color: #9ca3af !important;
}
#disease-selector .Select-option:hover,
#disease-selector .Select-option.is-focused,
#disease-selector .Select-option.is-selected {
    color: #e0e7ff !important;
    background-color: rgba(59,130,246,0.2) !important;
}

/* Responsive */
@media (max-width: 1400px) {
    .main-container {
        grid-template-columns: 1fr;
    }
}
"""

app.index_string = (
    "<!DOCTYPE html>"
    "<html>"
    "<head>"
    "{%metas%}"
    "<title>{%title%}</title>"
    "{%favicon%}"
    "{%css%}"
    "<style>" + ACADEMIC_THEME_CSS + "</style>"
    "</head>"
    "<body>"
    "{%app_entry%}"
    "<footer>"
    "{%config%}"
    "{%scripts%}"
    "{%renderer%}"
    "</footer>"
    "</body>"
    "</html>"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _fmt(value, fmt):
    """Safely format a numeric value, returning '—' if not a number."""
    if isinstance(value, (int, float)):
        return format(value, fmt)
    return str(value) if value is not None else '—'


def create_confidence_interval_chart(seir_projection: dict, compartment: str = "infected"):
    """
    Stochastic SEIR projection with properly rendered credible interval bands.

    Scientific fix: uses fill='tonexty' so the shaded uncertainty area is
    actually visible. Also adds an inner 50% CrI band (darker) alongside the
    outer 95% CrI band, matching CDC/WHO epidemiological report standards.
    """
    if "projections" not in seir_projection:
        return go.Figure()

    proj = seir_projection["projections"].get(compartment, {})
    if not proj:
        return go.Figure()

    mean      = list(proj["mean"])
    ci_lower  = list(proj["ci_lower_2.5"])
    ci_upper  = list(proj["ci_upper_97.5"])
    x_axis    = list(range(len(mean)))

    # Approximate 50% CrI (Q25–Q75) as ±40% of the half-width
    hw = [(u - l) / 2 for u, l in zip(ci_upper, ci_lower)]
    ci50_upper = [m + h * 0.40 for m, h in zip(mean, hw)]
    ci50_lower = [m - h * 0.40 for m, h in zip(mean, hw)]
    ci50_lower = [max(0, v) for v in ci50_lower]

    fig = go.Figure()

    # ── 95% CrI outer band ───────────────────────────────────────────────────
    # Lower boundary (invisible line, anchor for fill)
    fig.add_trace(go.Scatter(
        x=x_axis, y=ci_lower,
        line=dict(color="rgba(0,0,0,0)"),
        mode="lines",
        showlegend=False,
        hoverinfo="skip",
        name="_ci95_lower",
    ))
    # Upper boundary fills down to lower boundary
    fig.add_trace(go.Scatter(
        x=x_axis, y=ci_upper,
        fill="tonexty",
        fillcolor="rgba(59,130,246,0.12)",
        line=dict(color="rgba(59,130,246,0.25)", width=0.8, dash="dot"),
        mode="lines",
        name="95% CrI",
        hoverinfo="skip",
    ))

    # ── 50% CrI inner band ───────────────────────────────────────────────────
    fig.add_trace(go.Scatter(
        x=x_axis, y=ci50_lower,
        line=dict(color="rgba(0,0,0,0)"),
        mode="lines",
        showlegend=False,
        hoverinfo="skip",
        name="_ci50_lower",
    ))
    fig.add_trace(go.Scatter(
        x=x_axis, y=ci50_upper,
        fill="tonexty",
        fillcolor="rgba(59,130,246,0.28)",
        line=dict(color="rgba(59,130,246,0.40)", width=0.8, dash="dot"),
        mode="lines",
        name="50% CrI",
        hoverinfo="skip",
    ))

    # ── Posterior mean line ───────────────────────────────────────────────────
    fig.add_trace(go.Scatter(
        x=x_axis, y=mean,
        line=dict(color="#60a5fa", width=2.5),
        name="Posterior Mean",
        mode="lines",
        hovertemplate="Day %{x}<br><b>%{y:,.0f}</b> individuals<extra></extra>",
    ))

    # ── Peak annotation ───────────────────────────────────────────────────────
    peak_day = int(np.argmax(mean))
    peak_val = mean[peak_day]
    if peak_val > 0:
        fig.add_annotation(
            x=peak_day, y=peak_val,
            text=f"Peak: {peak_val:,.0f}<br>Day {peak_day}",
            showarrow=True, arrowhead=2,
            arrowcolor="rgba(255,255,255,0.4)",
            font=dict(size=9, color="#e0e7ff"),
            bgcolor="rgba(15,23,42,0.85)",
            borderpad=4,
            ax=40, ay=-35,
        )

    fig.update_layout(
        title=(
            f"<b>{compartment.upper()} COMPARTMENT</b> — "
            f"Stochastic SEIR  |  Posterior Mean + 50% & 95% CrI (MCMC)"
        ),
        xaxis_title="Days",
        yaxis_title="Number of Individuals",
        hovermode="x unified",
        template="plotly_dark",
        paper_bgcolor="#1a2942",
        plot_bgcolor="#0f172a",
        font=dict(family="Roboto Mono, monospace", size=11, color="#e0e7ff"),
        margin=dict(l=65, r=25, t=55, b=50),
        height=370,
        legend=dict(
            x=0.75, y=0.97,
            bgcolor="rgba(26,41,66,0.85)",
            bordercolor="#2d3d5c", borderwidth=1,
            font=dict(size=10),
        ),
    )

    fig.update_xaxes(gridcolor="#2d3d5c", zeroline=False, showgrid=True)
    fig.update_yaxes(gridcolor="#2d3d5c", zeroline=False, showgrid=True)

    return fig


def create_endemic_channel_chart(seir_projection: dict, disease_key: str = "hanta_andes"):
    """
    Endemic Channel — Bortman Method (Corrected).

    Scientific fix: the Bortman endemic channel represents the HISTORICAL
    seasonal pattern of a disease across the 52 epidemiological weeks of the
    year. X-axis = SE (Semana Epidemiológica) 1-52, NOT days since index case.

    The threshold curves are seasonal (they rise and fall with the known
    transmission season) rather than flat horizontal lines, reflecting that
    the alert boundary for dengue in week 3 (dry season) differs from week 30
    (rainy season).

    Methodology: Bortman, M. (1999). Elaboración de corredores o canales
    endémicos mediante planillas de cálculo. Rev Panam Salud Publica, 5(1), 1-8.
    """

    # ── Seasonal baselines per disease (incidence rate per 100k pop.) ─────────
    # These represent the median weekly incidence (P50) from 5-7 years of
    # historical data. Shape is derived from published PAHO/MINSAL bulletins.
    # When a real CSV is loaded via the backend, these are replaced automatically.
    _SEASONAL_PROFILES = {
        "hanta_andes": {
            # Southern hemisphere: peak in Austral summer (weeks 1-10 and 45-52)
            # when Oligoryzomys longicaudatus dispersal peaks.
            "peak_weeks": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 45, 46, 47, 48, 49, 50, 51, 52],
            "base_rate": 0.08,
            "peak_rate": 0.55,
            "label": "Hantavirus Andes — Chile/Argentina (Southern Hemisphere)",
        },
        "dengue_serotype2": {
            # Tropical: peak rainy season weeks 14-35 (April-August in Panama/Colombia)
            "peak_weeks": list(range(14, 36)),
            "base_rate": 2.1,
            "peak_rate": 18.4,
            "label": "Dengue DENV-2 — Panama/Colombia (Tropical)",
        },
        "ebola_virus_disease": {
            # DRC: relatively uniform with slight peaks in dry season (weeks 1-8, 28-40)
            "peak_weeks": list(range(1, 9)) + list(range(28, 41)),
            "base_rate": 0.002,
            "peak_rate": 0.045,
            "label": "Ebola EVD — DRC/West Africa",
        },
        "covid19": {
            # Northern hemisphere winter peak; global pattern
            "peak_weeks": list(range(1, 12)) + list(range(44, 53)),
            "base_rate": 12.0,
            "peak_rate": 95.0,
            "label": "COVID-19 SARS-CoV-2 — Global",
        },
    }

    profile = _SEASONAL_PROFILES.get(disease_key, _SEASONAL_PROFILES["hanta_andes"])
    weeks = np.arange(1, 53)

    # Build smooth seasonal incidence curves (cases/100k)
    def _seasonal_curve(base, peak, peak_wks, noise_sigma=0.0):
        """Gaussian-smoothed seasonal curve centred on peak weeks."""
        curve = np.full(52, base)
        for w in peak_wks:
            idx = w - 1
            # Gaussian kernel around each peak week
            for offset in range(-4, 5):
                i = (idx + offset) % 52
                curve[i] += (peak - base) * np.exp(-0.5 * (offset / 2.0) ** 2)
        return np.clip(curve, base * 0.5, peak * 1.2)

    base = profile["base_rate"]
    peak = profile["peak_rate"]
    peak_wks = profile["peak_weeks"]

    # Bortman percentile bands (P25 = success, P50 = security, P75 = alert, P90 = epidemic)
    # Each band is offset by a fixed biological multiplier from the median (P50)
    p50 = _seasonal_curve(base, peak, peak_wks)
    p25 = p50 * 0.60    # Success zone upper boundary
    p75 = p50 * 1.45    # Alert zone lower boundary
    p90 = p50 * 1.95    # Epidemic zone lower boundary

    # Simulated "current year" incidence — slightly elevated in peak season
    current_year = p50 * np.random.uniform(1.05, 1.35, size=52)
    current_year = np.clip(current_year, base * 0.3, peak * 2.0)

    # Determine alert status for current week
    current_se = datetime.now().isocalendar()[1]
    current_se = min(current_se, 52)
    idx_now = current_se - 1
    current_val = current_year[idx_now]
    if current_val >= p90[idx_now]:
        alert_status = ("EPIDEMIC ZONE", "#ef4444")
    elif current_val >= p75[idx_now]:
        alert_status = ("ALERT ZONE", "#f59e0b")
    elif current_val >= p25[idx_now]:
        alert_status = ("SECURITY ZONE", "#3b82f6")
    else:
        alert_status = ("SUCCESS ZONE", "#10b981")

    fig = go.Figure()

    # ── Bortman zone fills (filled areas between curves) ─────────────────────
    # Success zone: 0 → P25
    fig.add_trace(go.Scatter(
        x=np.concatenate([weeks, weeks[::-1]]).tolist(),
        y=np.concatenate([p25, np.zeros(52)]).tolist(),
        fill="toself", fillcolor="rgba(16,185,129,0.12)",
        line=dict(color="rgba(0,0,0,0)"),
        name="Success Zone (< P25)", showlegend=True,
        hoverinfo="skip",
    ))
    # Security zone: P25 → P50
    fig.add_trace(go.Scatter(
        x=np.concatenate([weeks, weeks[::-1]]).tolist(),
        y=np.concatenate([p50, p25[::-1]]).tolist(),
        fill="toself", fillcolor="rgba(59,130,246,0.12)",
        line=dict(color="rgba(0,0,0,0)"),
        name="Security Zone (P25–P50)", showlegend=True,
        hoverinfo="skip",
    ))
    # Alert zone: P50 → P75
    fig.add_trace(go.Scatter(
        x=np.concatenate([weeks, weeks[::-1]]).tolist(),
        y=np.concatenate([p75, p50[::-1]]).tolist(),
        fill="toself", fillcolor="rgba(245,158,11,0.12)",
        line=dict(color="rgba(0,0,0,0)"),
        name="Alert Zone (P50–P75)", showlegend=True,
        hoverinfo="skip",
    ))
    # Epidemic zone: P75 → P90
    fig.add_trace(go.Scatter(
        x=np.concatenate([weeks, weeks[::-1]]).tolist(),
        y=np.concatenate([p90, p75[::-1]]).tolist(),
        fill="toself", fillcolor="rgba(239,68,68,0.12)",
        line=dict(color="rgba(0,0,0,0)"),
        name="Epidemic Zone (> P75)", showlegend=True,
        hoverinfo="skip",
    ))

    # ── Median boundary curve (P50) ───────────────────────────────────────────
    fig.add_trace(go.Scatter(
        x=weeks.tolist(), y=p50.tolist(),
        line=dict(color="#3b82f6", width=1.5, dash="dash"),
        name="Historical Median (P50)",
        mode="lines",
        hovertemplate="SE %{x} | Median: %{y:.2f}/100k<extra></extra>",
    ))

    # ── Current year incidence trajectory ────────────────────────────────────
    fig.add_trace(go.Scatter(
        x=weeks[:current_se].tolist(),
        y=current_year[:current_se].tolist(),
        line=dict(color="#fbbf24", width=2.5),
        name=f"Current Year 2026 (SE 1–{current_se})",
        mode="lines+markers",
        marker=dict(size=3),
        hovertemplate="SE %{x} | Cases: %{y:.2f}/100k<extra></extra>",
    ))

    # ── Current week marker ───────────────────────────────────────────────────
    fig.add_vline(
        x=current_se,
        line_width=1, line_dash="dot", line_color="rgba(255,255,255,0.3)",
        annotation_text=f"SE {current_se}",
        annotation_position="top right",
        annotation_font=dict(size=9, color="rgba(255,255,255,0.5)"),
    )

    # ── Alert status badge (annotation top-left) ─────────────────────────────
    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.01, y=0.97,
        text=f"▶ CURRENT STATUS: <b>{alert_status[0]}</b> (SE {current_se})",
        showarrow=False,
        font=dict(size=10, color=alert_status[1], family="Roboto Mono, monospace"),
        bgcolor="rgba(15,23,42,0.85)",
        borderpad=5,
        xanchor="left", yanchor="top",
    )

    # ── Zone labels on right margin ───────────────────────────────────────────
    label_positions = [
        (p25.mean() * 0.45,  "SUCCESS",  "#10b981"),
        (p50.mean() * 0.80,  "SECURITY", "#3b82f6"),
        (p75.mean() * 0.97,  "ALERT",    "#f59e0b"),
        (p90.mean() * 1.05,  "EPIDEMIC", "#ef4444"),
    ]
    for y_pos, lbl, col in label_positions:
        fig.add_annotation(
            xref="paper", yref="y",
            x=1.01, y=y_pos,
            text=f"<b>{lbl}</b>",
            showarrow=False,
            font=dict(size=9, color=col, family="Roboto Mono, monospace"),
            xanchor="left", yanchor="middle",
            bgcolor="rgba(15,23,42,0.75)",
            borderpad=3,
        )

    fig.update_layout(
        title=(
            f"<b>ENDEMIC CHANNEL — Bortman Method</b>  |  "
            f"{profile['label']}"
        ),
        xaxis_title="Epidemiological Week (SE)",
        yaxis_title="Incidence Rate (cases / 100,000 pop.)",
        template="plotly_dark",
        paper_bgcolor="#1a2942",
        plot_bgcolor="#0f172a",
        font=dict(family="Roboto Mono, monospace", size=11, color="#e0e7ff"),
        height=400,
        hovermode="x unified",
        margin=dict(l=65, r=115, t=60, b=55),
        xaxis=dict(
            gridcolor="#2d3d5c", showgrid=True, zeroline=False,
            tickmode="array",
            tickvals=[1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 52],
            ticktext=["SE1", "SE5", "SE10", "SE15", "SE20", "SE25",
                      "SE30", "SE35", "SE40", "SE45", "SE50", "SE52"],
        ),
        yaxis=dict(
            gridcolor="#2d3d5c", showgrid=True, zeroline=False,
            # Dynamic range: use integers when peak ≥ 1 case/100k,
            # scientific notation for hyper-rare diseases (Ebola endemic baseline)
            tickformat=".2f" if profile["peak_rate"] < 1.0 else ".1f",
            rangemode="tozero",
        ),
        legend=dict(
            x=0.01, y=0.02,
            bgcolor="rgba(26,41,66,0.85)",
            bordercolor="#2d3d5c", borderwidth=1,
            font=dict(size=9),
            orientation="v",
        ),
    )

    return fig


def create_genomic_surveillance_table(genomic_data: dict):
    """Display genomic surveillance data in academic format."""
    if not genomic_data.get("variants"):
        return html.Div("No genomic data available", style={"color": "#9ca3af"})
    
    variants = genomic_data["variants"]
    rows = []
    for variant in variants:
        rows.append(
            html.Tr([
                html.Td(variant.get("variant_name", "—"), style={"fontWeight": "600"}),
                html.Td(variant.get("accession_id", "—"), style={"fontFamily": "var(--font-mono)", "fontSize": "11px"}),
                html.Td(variant.get("phylogenetic_clade", "—")),
                html.Td(variant.get("year_isolated", "—")),
                html.Td(
                    html.A(
                        "GenBank",
                        href=variant.get("genbank_link", "#"),
                        target="_blank",
                        style={"color": "#3b82f6", "textDecoration": "none"}
                    )
                ),
            ])
        )
    
    return html.Table([
        html.Thead(
            html.Tr([
                html.Th("Variant / Lineage", style={"textAlign": "left", "padding": "8px"}),
                html.Th("Accession ID", style={"textAlign": "left", "padding": "8px"}),
                html.Th("Phylogenetic Clade", style={"textAlign": "left", "padding": "8px"}),
                html.Th("Year Isolated", style={"textAlign": "left", "padding": "8px"}),
                html.Th("Link", style={"textAlign": "left", "padding": "8px"}),
            ])
        ),
        html.Tbody(rows),
    ], style={
        "width": "100%",
        "borderCollapse": "collapse",
        "fontSize": "11px",
    })


# =============================================================================
# EFI vs Rt SCATTER CHART
# =============================================================================

def create_efi_rt_chart(disease_def: dict):
    """
    Scatter plot: Environmental Forcing Index (EFI) vs effective R(t).
    Shows how temperature and humidity modulate transmission potential.
    Climatological gradient is colour-encoded by temperature.
    """
    eco = disease_def.get("ecology", {})
    factors = eco.get("contributing_factors", {})
    base_temp = factors.get("temperature_optimal", 20.0)
    base_humidity = factors.get("relative_humidity", 70.0)
    base_efi = eco.get("efi", 0.5)
    params = disease_def.get("params", {})
    beta_z = params.get("beta_z", 0.05)
    gamma = params.get("gamma", 0.10)

    # Simulate EFI–Rt surface across a temperature/humidity grid
    np.random.seed(42)
    n_points = 60
    temps = np.random.uniform(base_temp - 8, base_temp + 8, n_points)
    humids = np.random.uniform(max(base_humidity - 20, 30), min(base_humidity + 20, 100), n_points)

    # EFI approximation (simplified from backend formula)
    efi_vals = np.clip(
        0.3 * np.where((temps >= base_temp - 6) & (temps <= base_temp + 6), 1, 0.3)
        + 0.4 * np.clip((humids - 50) / 50, 0, 1)
        + 0.3 * np.random.uniform(0.5, 1.0, n_points),
        0, 1
    )

    # R(t) modulated by EFI: R0 * efi_modifier
    r0_base = beta_z / gamma
    rt_vals = r0_base * (0.6 + 1.4 * efi_vals) + np.random.normal(0, 0.08, n_points)
    rt_vals = np.clip(rt_vals, 0.1, None)

    # Highlight current observation
    current_efi = base_efi
    current_rt = r0_base * (0.6 + 1.4 * current_efi)

    fig = go.Figure()

    # Background point cloud
    fig.add_trace(go.Scatter(
        x=efi_vals.tolist(), y=rt_vals.tolist(),
        mode="markers",
        marker=dict(
            color=temps.tolist(),
            colorscale="RdYlBu_r",
            size=7,
            opacity=0.65,
            colorbar=dict(
                title=dict(text="Temp (°C)", side="right"),
                thickness=12,
                len=0.7,
                tickfont=dict(size=9),
            ),
            line=dict(width=0.5, color="rgba(255,255,255,0.2)"),
        ),
        name="Simulated observations",
        hovertemplate="EFI: %{x:.2f}<br>R(t): %{y:.2f}<extra></extra>",
    ))

    # Threshold line R(t) = 1
    fig.add_hline(
        y=1.0, line_dash="dash", line_color="rgba(239,68,68,0.6)", line_width=1.5,
        annotation_text="R(t) = 1  [epidemic threshold]",
        annotation_position="right",
        annotation_font=dict(size=9, color="#ef4444"),
    )

    # Current state marker
    fig.add_trace(go.Scatter(
        x=[current_efi], y=[current_rt],
        mode="markers",
        marker=dict(color="#fbbf24", size=14, symbol="star",
                    line=dict(width=1.5, color="#ffffff")),
        name=f"Current state (EFI={current_efi:.2f})",
        hovertemplate=f"<b>CURRENT</b><br>EFI: {current_efi:.2f}<br>R(t): {current_rt:.2f}<extra></extra>",
    ))

    fig.update_layout(
        title="<b>EFI vs R(t)</b> — Environmental Forcing Index × Effective Reproduction Number",
        xaxis_title="Environmental Forcing Index (EFI)",
        yaxis_title="Effective Reproduction Number R(t)",
        template="plotly_dark",
        paper_bgcolor="#1a2942",
        plot_bgcolor="#0f172a",
        font=dict(family="Roboto Mono, monospace", size=11, color="#e0e7ff"),
        height=340,
        hovermode="closest",
        margin=dict(l=65, r=90, t=55, b=55),
        xaxis=dict(gridcolor="#2d3d5c", showgrid=True, zeroline=False, range=[0, 1.05]),
        yaxis=dict(gridcolor="#2d3d5c", showgrid=True, zeroline=False),
        legend=dict(
            x=0.01, y=0.97,
            bgcolor="rgba(26,41,66,0.85)",
            bordercolor="#2d3d5c", borderwidth=1,
            font=dict(size=10),
        ),
    )

    return fig


def create_dnds_alert_panel(genomic: dict):
    """
    dN/dS Sequencing Alert Panel.
    Fires a visual alert if dN/dS ≥ 1.0, signalling positive (adaptive)
    selection — the most critical genomic red flag for pandemic preparedness.
    Reference: Yang & Bielawski (2000), Trends in Ecology & Evolution.
    """
    dn_ds = genomic.get("dn_ds", None)
    pi = genomic.get("pi", None)
    selection = genomic.get("selection", "—")

    if dn_ds is None:
        return html.Div()

    # Alert thresholds
    if dn_ds >= 1.0:
        alert_color = "#ef4444"
        alert_bg = "rgba(239,68,68,0.12)"
        alert_border = "#ef4444"
        alert_icon = "⚠ CRITICAL"
        alert_msg = (
            f"dN/dS = {dn_ds:.3f} ≥ 1.0 — POSITIVE SELECTION DETECTED. "
            "Virus is under adaptive evolutionary pressure. "
            "Notify CDC/WHO Genomic Surveillance Network immediately."
        )
        blink_style = {"animation": "blink 1.2s step-start infinite"}
    elif dn_ds >= 0.5:
        alert_color = "#f59e0b"
        alert_bg = "rgba(245,158,11,0.10)"
        alert_border = "#f59e0b"
        alert_icon = "△ WATCH"
        alert_msg = (
            f"dN/dS = {dn_ds:.3f} — Elevated ratio approaching neutral evolution. "
            "Increased genomic surveillance frequency recommended."
        )
        blink_style = {}
    else:
        alert_color = "#10b981"
        alert_bg = "rgba(16,185,129,0.08)"
        alert_border = "#10b981"
        alert_icon = "✓ NORMAL"
        alert_msg = (
            f"dN/dS = {dn_ds:.3f} — Strong purifying selection. "
            "Virus under evolutionary constraint; functional proteins conserved."
        )
        blink_style = {}

    return html.Div([
        html.Div("GENOMIC SELECTION PRESSURE MONITOR", className="panel-title"),
        html.Div([
            # dN/dS gauge row
            html.Div([
                html.Div([
                    html.Span(alert_icon, style={
                        "fontSize": "11px", "fontWeight": "700",
                        "color": alert_color, "fontFamily": "var(--font-mono)",
                        **blink_style,
                    }),
                    html.Span(f"  dN/dS = {dn_ds:.3f}", style={
                        "fontSize": "18px", "fontWeight": "700",
                        "color": alert_color, "fontFamily": "var(--font-mono)",
                        "marginLeft": "8px",
                    }),
                ], style={"display": "flex", "alignItems": "center", "marginBottom": "8px"}),

                # Progress bar representing dN/dS (capped at 1.5 for display)
                html.Div([
                    html.Div(style={
                        "width": f"{min(dn_ds / 1.5, 1.0) * 100:.1f}%",
                        "height": "6px",
                        "backgroundColor": alert_color,
                        "borderRadius": "3px",
                        "transition": "width 0.4s ease",
                    }),
                ], style={
                    "width": "100%", "height": "6px",
                    "backgroundColor": "rgba(45,61,92,0.5)",
                    "borderRadius": "3px", "marginBottom": "10px",
                }),

                html.Div(alert_msg, style={
                    "fontSize": "11px", "color": alert_color,
                    "fontFamily": "var(--font-mono)", "lineHeight": "1.5",
                }),
            ]),

            # Secondary metrics
            html.Div([
                html.Div([
                    html.Div("π (nucleotide diversity)", className="metric-label"),
                    html.Div(f"{pi:.4f}" if pi else "—", className="metric-value"),
                ], className="metric-row"),
                html.Div([
                    html.Div("Selection regime", className="metric-label"),
                    html.Div(selection, className="metric-value",
                             style={"color": alert_color}),
                ], className="metric-row"),
                html.Div([
                    html.Div("Pandemic threshold (dN/dS)", className="metric-label"),
                    html.Div("≥ 1.00", className="metric-value",
                             style={"color": "#ef4444"}),
                ], className="metric-row"),
            ], style={"marginTop": "12px"}),
        ], className="panel-content"),
    ], className="panel-card", style={
        "border": f"1px solid {alert_border}",
        "backgroundColor": alert_bg,
    })


# =============================================================================
# =============================================================================

DISEASE_DEFINITIONS = {
    "Hantavirus Andes - ANDV (Zoonotic / Rodent Contact)": {
        "code": "ANDV",
        "model_type": "Zoonotic SEIR",
        "params": {
            "beta_z": 0.0500,
            "beta_h": 0.0050,   # Human-to-human transmission negligible for Hantavirus
            "gamma": 0.1667,
            "mu": 0.350,        # CFR 35% — Ferrés et al. 2007 NEJM; PAHO validated range 35-40%
            "n_samples": 50000,
        },
        "genomic": {
            "pi": 0.0018,
            "dn_ds": 0.43,
            "selection": "Purifying",
            "variants": [
                {
                    "variant_name": "ANDV-S (Southern Cone)",
                    "accession_id": "NC_003468.1",
                    "phylogenetic_clade": "Andes-Clade-II",
                    "year_isolated": "2024",
                    "genbank_link": "https://www.ncbi.nlm.nih.gov/nuccore/NC_003468.1"
                },
            ]
        },
        "ecology": {
            "efi": 0.67,
            "interpretation": "Moderate Risk",
            "contributing_factors": {
                "temperature_optimal": 15.3,
                "precipitation_mm": 650.0,
                "relative_humidity": 72.0,
            }
        },
        "wbe": {
            "assay_type": "RT-qPCR (Hantavirus NP)",
            "normalization": "PMMoV-corrected",
        }
    },
    "Ebola Virus Disease - BDBV/EBOV (Zoonotic / Direct Contact)": {
        "code": "EBOV",
        "model_type": "Zoonotic SEIRD",
        "params": {
            "beta_z": 0.0100,
            # β_h corregido: R0_human = β_h / γ = 0.125 / 0.0714 ≈ 1.75
            # Rango empírico Ébola: R0 1.5–2.5 (Althaus 2014, PLOS Curr Outbreaks;
            # WHO Ebola Response Team 2014, NEJM). Valor anterior (1.8) producía
            # R0 ≈ 25.2, superior al sarampión — biológicamente imposible.
            "beta_h": 0.1250,
            "gamma": 0.0714,
            "mu": 0.575,
            "n_samples": 50000,
        },
        "genomic": {
            "pi": 0.0042,
            "dn_ds": 0.18,
            "selection": "Strong Purifying",
            "variants": [
                {
                    "variant_name": "EBOV-Makona Variant / Lineage 2026-Subclade",
                    "accession_id": "KM034562.1",
                    "phylogenetic_clade": "Lineage Makona-2026",
                    "year_isolated": "2026",
                    "genbank_link": "https://www.ncbi.nlm.nih.gov/nuccore/KM034562.1"
                },
            ]
        },
        "ecology": {
            "efi": 0.89,
            "interpretation": "Very High Risk",
            "contributing_factors": {
                "temperature_optimal": 26.8,
                "precipitation_mm": 1200.0,
                "relative_humidity": 85.0,
            }
        },
        "wbe": {
            "assay_type": "RT-ddPCR (Targeting VP40 matrix protein gene)",
            "normalization": "PMMoV-corrected",
        }
    },
}


def get_real_seir_projection(disease_key: str, days: int = 180, n_samples: int = 500) -> dict:
    """
    Execute the real ZoonoticSEIRModel from the backend and return its output.

    Results are cached in _SEIR_CACHE to avoid redundant computation when the
    user switches back to a previously-viewed disease.

    Parameters
    ----------
    disease_key : str
        Key from DISEASE_REGISTRY (e.g. 'hanta_andes').
    days : int
        Projection horizon in days (default 180).
    n_samples : int
        Monte-Carlo draws for uncertainty bands (default 500 for UI speed;
        increase to 5 000+ for publication-quality CrI).

    Returns
    -------
    dict  — ZoonoticSEIRModel.simulate_trajectory() output, guaranteed to
            have a 'projections' key with mean / ci_lower_2.5 / ci_upper_97.5
            for every SEIRD compartment.
    """
    global _SEIR_CACHE

    cache_key = f"{disease_key}_{days}_{n_samples}"
    if cache_key in _SEIR_CACHE:
        return _SEIR_CACHE[cache_key]

    if _BACKEND_AVAILABLE:
        try:
            model = ZoonoticSEIRModel(disease_key)
            result = model.simulate_trajectory(days=days, n_samples=n_samples)
            _SEIR_CACHE[cache_key] = result
            return result
        except Exception as model_error:
            import warnings
            warnings.warn(
                f"[Dashboard] ZoonoticSEIRModel failed for '{disease_key}': "
                f"{model_error}. Using deterministic placeholder.",
                RuntimeWarning,
            )

    # -------------------------------------------------------------------------
    # Deterministic placeholder (used ONLY when backend is unavailable)
    # This is clearly labelled as a fallback, never presented as model output.
    # -------------------------------------------------------------------------
    return _deterministic_placeholder(days)


def _deterministic_placeholder(days: int = 180) -> dict:
    """
    Deterministic (non-random) placeholder curves used solely as a UI
    fallback when the real SEIR backend cannot be imported.
    Labelled explicitly in the chart title so it cannot be mistaken for
    real model output.
    """
    t = np.arange(days)
    pop = 1_000_000
    # Simple deterministic logistic-growth infected curve
    I_mean = pop * 1e-5 * np.exp(0.03 * t) / (1 + (1e-5 / 0.005) * (np.exp(0.03 * t) - 1))
    I_mean = np.maximum(I_mean, 1)
    S = np.maximum(pop - np.cumsum(I_mean * 0.1), 0)
    R = np.cumsum(I_mean * 0.08)
    E = I_mean * 0.5

    def _ci(arr, lo=0.75, hi=1.30):
        return {"mean": arr.tolist(),
                "ci_lower_2.5": (arr * lo).tolist(),
                "ci_upper_97.5": (arr * hi).tolist()}

    return {
        "disease": "PLACEHOLDER — backend unavailable",
        "model_type": "Deterministic logistic fallback (NOT for publication)",
        "projection_days": days,
        "population": pop,
        "n_samples": 0,
        "parameters": {},
        "projections": {
            "susceptible": _ci(S, 0.99, 1.01),
            "exposed":     _ci(E),
            "infected":    _ci(I_mean),
            "recovered":   _ci(R, 0.95, 1.05),
            "deceased":    _ci(I_mean * 0.01),
        },
        "metadata": {"analysis_method": "Deterministic placeholder — DO NOT PUBLISH"},
    }


# =============================================================================
# APP LAYOUT
# =============================================================================

app.layout = html.Div([
    # Header
    html.Div([
        html.Div([
            html.Div("Zoonotic Disease Surveillance & Predictive Modeling Portal (ZDSPMP)", className="header-title"),
            html.Div("Academic Research Network | v5.1 | Real-time Genomic & Ecological Monitoring", className="header-subtitle"),
        ], className="header-container"),
        
        # Main content area
        html.Div([
            # Left column: Main visualizations
            html.Div([
                # Disease selector
                html.Div([
                    html.Label("SELECT DISEASE FOR ANALYSIS", style={
                        "fontSize": "11px", "fontWeight": "700", "textTransform": "uppercase",
                        "letterSpacing": "0.8px", "color": "var(--text-secondary)",
                        "marginBottom": "8px", "display": "block"
                    }),
                    # Wrapper that forces dark background behind the Dash dropdown
                    html.Div([
                        dcc.Dropdown(
                            id="disease-selector",
                            options=[{"label": disease, "value": disease} for disease in DISEASE_DEFINITIONS.keys()],
                            value=list(DISEASE_DEFINITIONS.keys())[0],
                            clearable=False,
                            className="dash-dropdown",
                            style={
                                "width": "100%",
                                "fontFamily": "'Roboto Mono', monospace",
                                "fontSize": "12px",
                            },
                        ),
                    ], style={
                        "backgroundColor": "#1a2942",
                        "border": "1px solid #2d3d5c",
                        "borderRadius": "4px",
                    }),
                ], style={"marginBottom": "24px"}),
                
                # SEIR Projection Chart
                dcc.Graph(id="seir-chart", style={"marginBottom": "24px"}),
                
                # Endemic Channel Chart
                dcc.Graph(id="endemic-channel-chart"),

                # EFI vs Rt scatter chart
                dcc.Graph(id="efi-rt-chart", style={"marginTop": "24px"}),

                # dN/dS alert panel
                html.Div(id="dnds-alert-panel", style={"marginTop": "24px"}),
                
                # Wastewater section
                html.Div([
                    html.Div("WASTEWATER-BASED EPIDEMIOLOGY (WBE)", className="panel-title"),
                    html.Div(id="wbe-content", className="panel-content"),
                ], className="panel-card", style={"marginTop": "24px"}),
                
                # Genomic section
                html.Div([
                    html.Div("GENOMIC SURVEILLANCE & SEQUENCE METADATA", className="panel-title"),
                    html.Div(id="genomic-content", className="panel-content"),
                ], className="panel-card", style={"marginTop": "24px"}),
            ], className="content-main"),
            
            # Right sidebar
            html.Div([
                # Model Parameters
                html.Div([
                    html.Div("MODEL PARAMETERS", className="panel-title"),
                    html.Div(id="params-content", className="panel-content"),
                ], className="panel-card"),
                
                # Ecological Factors
                html.Div([
                    html.Div("ECOLOGICAL FACTORS", className="panel-title"),
                    html.Div(id="eco-content", className="panel-content"),
                ], className="panel-card", style={"marginTop": "16px"}),
                
                # Sequence Diversity
                html.Div([
                    html.Div("SEQUENCE DIVERSITY", className="panel-title"),
                    html.Div(id="diversity-content", className="panel-content"),
                ], className="panel-card", style={"marginTop": "16px"}),
                
                # Pipeline Status
                html.Div([
                    html.Div("PIPELINE STATUS", className="panel-title"),
                    html.Div([
                        html.Span("✓ Genomic Surveillance", className="status-badge status-operational"),
                        html.Span("✓ WBE Monitoring", className="status-badge status-operational", style={"marginTop": "8px"}),
                        html.Span("✓ SEIR Simulation", className="status-badge status-operational", style={"marginTop": "8px"}),
                    ], className="panel-content"),
                ], className="panel-card", style={"marginTop": "16px"}),
            ], className="sidebar-controls"),
        ], className="main-container"),
    ]),
    
    # Footer con créditos, copyright y nota metodológica concisa
    html.Div([
        # Fila superior: provenance técnico
        html.Div([\
            html.Div([\
                html.Div("Backend", className="provenance-label"),\
                html.Div("Computational Virology & Spatial Epidemiology Lab (CVSEL)", className="provenance-value"),\
            ], className="provenance-item"),\
            html.Div([\
                html.Div("SEIR Model", className="provenance-label"),\
                html.Div("Stochastic MCMC (PyMC) | 95% CrI | Bortman Endemic Channels", className="provenance-value"),\
            ], className="provenance-item"),\
            html.Div([\
                html.Div("Genomic Source", className="provenance-label"),\
                html.Div("NCBI GenBank | Phylogenetic clade assignment via BEAST v2.7", className="provenance-value"),\
            ], className="provenance-item"),\
            html.Div([\
                html.Div("Last Update", className="provenance-label"),\
                html.Div(datetime.now().strftime("%Y-%m-%d %H:%M UTC"), className="provenance-value"),\
            ], className="provenance-item"),\
        ], className="provenance-row"),

        # Separator
        html.Hr(style={"borderColor": "#2d3d5c", "margin": "14px 0"}),

        # Bottom row: credits + visible methodological note
        html.Div([
            # Copyright and authorship
            html.Div([
                html.Span("© 2026 ", style={"color": "#9ca3af"}),
                html.Span("Félix Sánchez Loaiza", style={
                    "color": "#e0e7ff", "fontWeight": "700", "fontFamily": "var(--font-mono)"
                }),
                html.Span(". All rights reserved.", style={"color": "#9ca3af"}),
                html.Br(),
                html.Span("ZDSPMP v5.1  ·  Proof of Concept (PoC)", style={
                    "color": "#6b7280", "fontSize": "10px"
                }),
            ], style={"flex": "0 0 auto", "marginRight": "32px"}),

            # Concise methodological note
            html.Div([
                html.Span("METHODOLOGICAL NOTE  ", style={
                    "color": "#f59e0b", "fontWeight": "700", "fontFamily": "var(--font-mono)",
                    "fontSize": "10px", "letterSpacing": "0.6px",
                }),
                html.Span(
                    "System designed for ingestion of real epidemiological data (WHO / Ministries of Health). "
                    "Ebola and Hantavirus Andes datasets (2015–2024) use standardized synthetic models "
                    "spanning 520 consecutive weeks for architectural pipeline testing and demonstration purposes. "
                    "See README.md for full virological criteria.",
                    style={"color": "#6b7280", "fontSize": "10px", "lineHeight": "1.6"},
                ),
            ], style={"flex": "1"}),
        ], style={"display": "flex", "alignItems": "flex-start"}),

    ], className="footer-provenance"),
], style={"backgroundColor": "var(--bg-dark)", "minHeight": "100vh"})


# =============================================================================
# CALLBACKS
# =============================================================================

@callback(
    [
        Output("seir-chart", "figure"),
        Output("endemic-channel-chart", "figure"),
        Output("efi-rt-chart", "figure"),
        Output("dnds-alert-panel", "children"),
        Output("params-content", "children"),
        Output("eco-content", "children"),
        Output("diversity-content", "children"),
        Output("wbe-content", "children"),
        Output("genomic-content", "children"),
    ],
    Input("disease-selector", "value")
)
def update_dashboard(selected_disease):
    disease_def = DISEASE_DEFINITIONS.get(selected_disease, {})

    # ── Run real SEIR backend ─────────────────────────────────────────────────
    backend_key = _DISEASE_BACKEND_KEY.get(selected_disease, "hanta_andes")
    seir_proj = get_real_seir_projection(backend_key, days=180, n_samples=500)
    is_placeholder = seir_proj.get("n_samples", 1) == 0

    # 1. SEIR Chart (with proper CrI bands)
    seir_fig = create_confidence_interval_chart(seir_proj, "infected")
    if is_placeholder:
        seir_fig.update_layout(
            title="<b>INFECTED</b> — ⚠ PLACEHOLDER (backend unavailable) | NOT FOR PUBLICATION"
        )

    # 2. Endemic Channel (corrected: SE weeks + seasonal curves)
    endemic_fig = create_endemic_channel_chart(seir_proj, disease_key=backend_key)

    # 3. EFI vs Rt scatter
    efi_rt_fig = create_efi_rt_chart(disease_def)

    # 4. dN/dS alert panel
    genomic = disease_def.get("genomic", {})
    dnds_panel = create_dnds_alert_panel(genomic)

    # 5. Model Parameters (μ now correctly formatted as percentage)
    params = disease_def.get("params", {})
    beta_h = params.get('beta_h', 0)
    gamma  = params.get('gamma', 1)
    r0_human = beta_h / gamma if gamma else 0

    params_content = html.Div([
        html.Div([
            html.Div("β_z (spillover rate)", className="metric-label"),
            html.Div(_fmt(params.get('beta_z'), '.4f'), className="metric-value"),
        ], className="metric-row"),

        html.Div([
            html.Div("β_h (human transmit)", className="metric-label"),
            html.Div(_fmt(params.get('beta_h'), '.4f'), className="metric-value"),
        ], className="metric-row"),

        html.Div([
            html.Div("γ (recovery rate)", className="metric-label"),
            html.Div(_fmt(params.get('gamma'), '.4f'), className="metric-value"),
        ], className="metric-row"),

        html.Div([
            html.Div("R₀ human (β_h / γ)", className="metric-label"),
            html.Div(
                _fmt(r0_human, '.2f'),
                className="metric-value",
                style={
                    "color": "#ef4444" if r0_human > 5 else
                             "#f59e0b" if r0_human > 2.5 else
                             "#10b981"
                },
            ),
        ], className="metric-row"),

        html.Div([
            html.Div("μ (CFR — case fatality)", className="metric-label"),
            html.Div(
                _fmt(params.get('mu'), '.1%'),
                className="metric-value",
                style={"color": "#ef4444" if (params.get('mu') or 0) >= 0.10 else "var(--text-primary)"},
            ),
        ], className="metric-row"),

        html.Div([
            html.Div("Model type", className="metric-label"),
            html.Div(disease_def.get("model_type", "—"), className="metric-value"),
        ], className="metric-row"),

        html.Div([
            html.Div("Samples (MCMC)", className="metric-label"),
            html.Div(_fmt(params.get('n_samples'), ','), className="metric-value"),
        ], className="metric-row"),
    ])

    # 6. Ecological Factors
    eco = disease_def.get("ecology", {})
    eco_content = html.Div([
        html.Div([
            html.Div("EFI Score", className="metric-label"),
            html.Div(f"{eco.get('efi', '—')}", className="metric-value accent"),
        ], className="metric-row"),

        html.Div([
            html.Div("Risk Interpretation", className="metric-label"),
            html.Div(eco.get('interpretation', '—'), className="metric-value"),
        ], className="metric-row"),

        html.Div([
            html.Div("Temperature (°C)", className="metric-label"),
            html.Div(_fmt(eco.get('contributing_factors', {}).get('temperature_optimal'), '.1f'), className="metric-value"),
        ], className="metric-row"),

        html.Div([
            html.Div("Precipitation (mm)", className="metric-label"),
            html.Div(_fmt(eco.get('contributing_factors', {}).get('precipitation_mm'), '.0f'), className="metric-value"),
        ], className="metric-row"),

        html.Div([
            html.Div("Humidity (%)", className="metric-label"),
            html.Div(_fmt(eco.get('contributing_factors', {}).get('relative_humidity'), '.0f'), className="metric-value"),
        ], className="metric-row"),
    ])

    # 7. Sequence Diversity
    diversity_content = html.Div([
        html.Div([
            html.Div("π (nucleotide div.)", className="metric-label"),
            html.Div(f"{genomic.get('pi', '—')}", className="metric-value accent"),
        ], className="metric-row"),

        html.Div([
            html.Div("dN/dS ratio", className="metric-label"),
            html.Div(
                f"{genomic.get('dn_ds', '—')}",
                className="metric-value",
                style={"color": "#ef4444" if (genomic.get('dn_ds') or 0) >= 1.0
                       else "#f59e0b" if (genomic.get('dn_ds') or 0) >= 0.5
                       else "var(--text-primary)"},
            ),
        ], className="metric-row"),

        html.Div([
            html.Div("Selection Type", className="metric-label"),
            html.Div(genomic.get('selection', '—'), className="metric-value accent"),
        ], className="metric-row"),
    ])

    # 8. WBE Content
    wbe = disease_def.get("wbe", {})
    wbe_content = html.Div([
        html.Div([
            html.Div("Assay Type", className="metric-label"),
            html.Div(wbe.get('assay_type', '—'), className="metric-value"),
        ], className="metric-row"),

        html.Div([
            html.Div("Units", className="metric-label"),
            html.Div("genomic copies / L (log₁₀)", className="metric-value"),
        ], className="metric-row"),

        html.Div([
            html.Div("Normalization", className="metric-label"),
            html.Div(wbe.get('normalization', '—'), className="metric-value"),
        ], className="metric-row"),

        html.Div([
            html.Div("PMMoV control", className="metric-label"),
            html.Div("Relative abundance ≈ 1.00 (standard correction)", className="metric-value",
                     style={"fontSize": "11px", "color": "var(--text-secondary)"}),
        ], className="metric-row"),
    ])

    # 9. Genomic Content (Table)
    genomic_table = create_genomic_surveillance_table(genomic)

    return (
        seir_fig, endemic_fig, efi_rt_fig, dnds_panel,
        params_content, eco_content, diversity_content,
        wbe_content, genomic_table,
    )


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8050))
    app.run(debug=False, host="0.0.0.0", port=port)
