import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import streamlit as st

from utils.state import init_state
from utils.ui import setup_sidebar
from core.plots import fig_to_png
from core.registry import add_plot

# Draw every point (no line simplification)
mpl.rcParams["path.simplify"] = False

st.set_page_config(page_title="MINT LLM", page_icon=":herb:", layout="wide")
init_state()
setup_sidebar()

st.title("🧪 Simulation Assessment (Logs)")
st.caption(
    "Upload MD log files; we’ll parse any `Step`-based tables and auto-plot "
    "every numeric signal vs the first numeric x-axis (often Step or Time)."
)

# ---------- Session storage for logs & de-dup of Results ----------
if "log_state" not in st.session_state:
    st.session_state.log_state = {
        "engine": "",
        "files": [],                 # list[(name, bytes)]
        "tables": {},                # name -> DataFrame
        "plots": {},                 # name -> { ycol -> png_bytes }
        "added_to_results": set(),   # set of "name|ycol"
    }

engine = st.session_state.get("engine", "")
if not engine:
    st.info("Pick an engine in the left sidebar to proceed.")
    st.stop()

# ---------- Page-level display options ----------
col_opt1, col_opt2 = st.columns(2)
with col_opt1:
    smoothing = st.checkbox("Show moving average (for readability)", value=False)
with col_opt2:
    window = st.number_input(
        "Smoothing window (points)",
        min_value=3, value=11, step=2,
        help="Used only if smoothing is enabled; odd number recommended."
    )

# ---------- File uploader (multiple) ----------
new_files = st.file_uploader(
    "Upload one or more log files",
    type=["log", "lammps", "mdout", "mdinfo", "txt"],
    accept_multiple_files=True
)

# ---------- Parser ----------
def parse_all_step_tables(text: str) -> pd.DataFrame:
    """
    Collect all whitespace-delimited tables whose header line starts with 'Step'
    and whose rows start with digits. Concatenate them, coerce numeric, sort,
    and drop duplicates.
    """
    lines = [ln.rstrip("\n") for ln in text.splitlines()]
    dfs = []
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln.startswith("Step "):  # header found
            headers = ln.split()
            i += 1
            rows = []
            while i < len(lines):
                row = lines[i].strip()
                if not row or not row[0].isdigit():
                    break
                parts = row.split()
                if len(parts) == len(headers):
                    rows.append(parts)
                else:
                    break
                i += 1
            if rows:
                df = pd.DataFrame(rows, columns=headers)
                for c in df.columns:
                    df[c] = pd.to_numeric(df[c], errors="ignore")
                dfs.append(df)
        else:
            i += 1

    if not dfs:
        return pd.DataFrame()
    out = pd.concat(dfs, ignore_index=True)

    # Sort by Step then Time if available; otherwise by the first numeric column
    sort_cols = []
    if "Step" in out.columns:
        sort_cols.append("Step")
    if "Time" in out.columns:
        sort_cols.append("Time")
    if not sort_cols:
        num_cols_tmp = [c for c in out.columns if pd.api.types.is_numeric_dtype(out[c])]
        if num_cols_tmp:
            sort_cols = [num_cols_tmp[0]]
    if sort_cols:
        out = out.sort_values(sort_cols).reset_index(drop=True)

    return out.drop_duplicates()

# ---------- Persistence logic ----------
if st.session_state.log_state["files"] and not new_files:
    st.success("Using previously uploaded log files (persisted).")
elif new_files:
    # Replace state when new uploads provided
    st.session_state.log_state = {
        "engine": engine,
        "files": [(f.name, f.getbuffer().tobytes()) for f in new_files],
        "tables": {},
        "plots": {},
        "added_to_results": set(),
    }

# ---------- Process & plot ----------
if st.session_state.log_state["files"]:
    st.markdown("### Current files")
    for (fname, _) in st.session_state.log_state["files"]:
        st.write(f"- `{fname}`")

    # Parse and generate plots for any files not yet processed
    for (fname, blob) in st.session_state.log_state["files"]:
        if fname in st.session_state.log_state["tables"]:
            continue

        text = blob.decode(errors="ignore")
        df = parse_all_step_tables(text)

        st.session_state.log_state["tables"][fname] = df
        st.session_state.log_state["plots"][fname] = {}

        if df.empty:
            continue

        # Determine numeric columns once per file
        num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        if len(num_cols) < 2:
            continue
        xcol = num_cols[0]

        # One decimation slider per file (unique key)
        decimate = st.slider(
            f"Plot every Nth point for `{fname}`",
            min_value=1, max_value=50, value=1, step=1,
            help="1 = plot all points",
            key=f"decimate_{fname}"
        )

        for ycol in num_cols[1:]:
            # Raw arrays
            x = pd.to_numeric(df[xcol], errors="coerce").values
            y = pd.to_numeric(df[ycol], errors="coerce").values

            # Start from raw; optional smoothing for display only
            y_plot = y.copy()
            if smoothing and window > 1 and window % 2 == 1 and len(y_plot) >= window:
                y_plot = pd.Series(y_plot).rolling(window, center=True, min_periods=1).mean().values

            # Decimate plot (not CSV)
            x_plot = x[::decimate]
            y_plot = y_plot[::decimate]

            # Plot to look like RMSD: thin line, grid, every point (no markers needed)
            fig, ax = plt.subplots()
            ax.plot(x_plot, y_plot, linestyle="-", linewidth=1.0, antialiased=True)
            ax.set_xlabel(xcol); ax.set_ylabel(ycol)
            ax.set_title(f"{ycol} vs {xcol}")
            ax.grid(True, alpha=0.25)
            png = fig_to_png(fig)

            # Persist for revisit
            st.session_state.log_state["plots"][fname][ycol] = png

            # CSV of raw data (non-decimated, non-smoothed)
            csv_buf = io.BytesIO()
            pd.DataFrame({xcol: x, ycol: y}).to_csv(csv_buf, index=False)
            csv_bytes = csv_buf.getvalue()

            # Push to Results once per (file,ycol)
            key = f"{fname}|{ycol}"
            if key not in st.session_state.log_state["added_to_results"]:
                add_plot(
                    st.session_state,
                    kind="log_timeseries",
                    title=f"{ycol} vs {xcol} [{fname}]",
                    params={"file": fname, "x": xcol, "y": ycol, "engine": engine,
                            "decimate": decimate, "smoothed": bool(smoothing)},
                    png_bytes=png,
                    csv_bytes=csv_bytes,
                    tags=[engine, "assessment"],
                )
                st.session_state.log_state["added_to_results"].add(key)

    # ---------- Show parsed tables + persisted plots ----------
    st.markdown("### Parsed results")
    for fname, df in st.session_state.log_state["tables"].items():
        st.subheader(fname)
        if df.empty:
            st.warning("No Step-based table found.")
            continue

        st.dataframe(df.head())

        plot_dict = st.session_state.log_state["plots"].get(fname, {})
        if plot_dict:
            cols = st.columns(2)
            i = 0
            for ycol, png in plot_dict.items():
                with cols[i % 2]:
                    st.image(png, caption=f"{ycol} vs {df.columns[0]} (auto-plot)")
                i += 1
        else:
            st.info("No numeric columns to plot.")

    # Controls
    if st.button("Refresh / Replace logs"):
        st.session_state.log_state = {
            "engine": engine,
            "files": [],
            "tables": {},
            "plots": {},
            "added_to_results": set(),
        }
        st.rerun()
else:
    st.info("Upload one or more log files to begin.")
