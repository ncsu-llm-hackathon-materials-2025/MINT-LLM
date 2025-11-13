import streamlit as st
from utils.state import init_state
from utils.ui import setup_sidebar
from utils.io_helpers import save_uploaded_file
from core.loaders import load_universe
import os

st.set_page_config(page_title="MINT LLM", page_icon=":herb:", layout="wide")
init_state()
setup_sidebar()

st.title("📦 Trajectory Upload")
st.caption("Upload topology + trajectory files for LAMMPS, GROMACS, or AMBER.")

# --- USE THE GLOBAL ENGINE FROM SIDEBAR ---
engine = st.session_state.get("engine", "")
if not engine:
    st.info("Pick an engine in the left sidebar to proceed.")
    st.stop()

eng_key = engine.lower()
if "uploads" not in st.session_state:
    st.session_state["uploads"] = {}
if eng_key not in st.session_state["uploads"]:
    st.session_state["uploads"][eng_key] = {}

uploaded = st.session_state["uploads"][eng_key]  # bind to session

# --- uploader UI ---
if engine == "LAMMPS":
    col1, col2 = st.columns(2)
    with col1:
        dataf = st.file_uploader("data.lammps", type=["lammps","data"])
        if dataf:
            uploaded["data.lammps"] = save_uploaded_file(dataf)
    with col2:
        dumps = st.file_uploader("dump.lammpstrj (1+)", type=["lammpstrj","dump"], accept_multiple_files=True)
        if dumps:
            uploaded["dump.lammpstrj"] = [save_uploaded_file(d) for d in dumps]

elif engine == "GROMACS":
    col1, col2 = st.columns(2)
    with col1:
        top = st.file_uploader("Topology (topol.tpr or .gro)", type=["tpr","gro"])
        if top:
            saved = save_uploaded_file(top)
            if top.name.lower().endswith(".tpr"):
                uploaded["topol.tpr"] = saved
                uploaded.pop("topol.gro", None)
            else:
                uploaded["topol.gro"] = saved
                uploaded.pop("topol.tpr", None)
    with col2:
        xtcs = st.file_uploader("traj (xtc/trr) 1+", type=["xtc","trr"], accept_multiple_files=True)
        if xtcs:
            uploaded["traj"] = [save_uploaded_file(x) for x in xtcs]

else:  # AMBER
    col1, col2 = st.columns(2)
    with col1:
        prmtop = st.file_uploader("prmtop", type=["prmtop","top","parm7"])
        if prmtop:
            uploaded["prmtop"] = save_uploaded_file(prmtop)
    with col2:
        ncs = st.file_uploader("prod.nc (1+)", type=["nc","netcdf"], accept_multiple_files=True)
        if ncs:
            uploaded["prod"] = [save_uploaded_file(n) for n in ncs]

# persist uploads
st.session_state["uploads"][eng_key] = uploaded

# show what's saved
if uploaded:
    st.markdown("**Current files saved for this engine:**")
    for k, v in uploaded.items():
        if isinstance(v, list):
            for i, p in enumerate(v, 1):
                st.write(f"- {k}[{i}]: `{os.path.basename(p)}`")
        else:
            st.write(f"- {k}: `{os.path.basename(v)}`")
    if st.button("Clear files for this engine"):
        st.session_state["uploads"][eng_key] = {}
        st.rerun()

# build Universe (DO NOT write to st.session_state['engine'] here)
if st.button("Create/Refresh Universe"):
    try:
        u = load_universe(engine, uploaded)
        st.session_state["universe_ref"] = u
        st.success(f"Universe created: {u.atoms.n_atoms} atoms, {u.trajectory.n_frames} frames")
    except Exception as e:
        st.error(f"Failed to create Universe: {e}")
