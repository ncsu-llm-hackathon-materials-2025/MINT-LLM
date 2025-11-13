import streamlit as st
import tempfile, os
import MDAnalysis as mda
from utils.state import init_state
from utils.ui import setup_sidebar

st.set_page_config(page_title="MINT LLM", page_icon=":herb:", layout="wide")
init_state()
setup_sidebar()
st.title("🔁 Trajectory Conversion")
st.caption("Convert to exchange formats (XTC/TRR/DCD/LAMMPS dump) and optional PDB/GRO snapshot.")

if not st.session_state.get("universe_ref"):
    st.warning("Upload a trajectory and create a Universe in Trajectory Upload first.")
else:
    u = st.session_state["universe_ref"]
    sel = st.text_input("Atom selection (MDAnalysis syntax)", value="all")
    stride = st.number_input("Frame stride", min_value=1, value=10, step=1)
    traj_ext = st.selectbox("Trajectory format", [".xtc",".trr",".dcd",".lammpstrj"], index=0)
    write_top = st.checkbox("Also write a single-frame topology snapshot")
    top_ext = st.selectbox("Topology format", [".pdb",".gro"], index=0) if write_top else None

    if st.button("Convert"):
        ag = u.select_atoms(sel)
        with tempfile.TemporaryDirectory() as td:
            traj_path = os.path.join(td, f"converted{traj_ext}")
            if traj_ext == ".xtc":
                with mda.coordinates.XTC.XTCWriter(traj_path, n_atoms=ag.n_atoms) as W:
                    for ts in u.trajectory[::stride]:
                        W.write(ag)
            elif traj_ext == ".trr":
                with mda.coordinates.TRR.TRRWriter(traj_path, n_atoms=ag.n_atoms) as W:
                    for ts in u.trajectory[::stride]:
                        W.write(ag)
            elif traj_ext == ".dcd":
                with mda.coordinates.DCD.DCDWriter(traj_path, n_atoms=ag.n_atoms) as W:
                    for ts in u.trajectory[::stride]:
                        W.write(ag)
            else:
                with mda.coordinates.LAMMPSDUMP.DUMPWriter(traj_path, n_atoms=ag.n_atoms) as W:
                    for ts in u.trajectory[::stride]:
                        W.write(ag)

            with open(traj_path, "rb") as f:
                st.download_button("Download trajectory", f.read(), file_name=f"converted{traj_ext}")

            if write_top:
                top_path = os.path.join(td, f"snapshot{top_ext}")
                ag.write(top_path)
                with open(top_path, "rb") as f:
                    st.download_button("Download topology snapshot", f.read(), file_name=f"snapshot{top_ext}")
