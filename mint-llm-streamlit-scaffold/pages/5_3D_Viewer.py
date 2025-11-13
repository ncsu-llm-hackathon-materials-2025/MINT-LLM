import streamlit as st
from utils.state import init_state
from utils.ui import setup_sidebar
import py3Dmol
import tempfile
import os

st.set_page_config(page_title="MINT LLM", page_icon=":herb:", layout="wide")
init_state()
setup_sidebar()

st.title("🧬 3D Viewer (Trajectory Animation)")
st.caption("Preview a subsampled trajectory with py3Dmol. Choose a selection and how many frames to animate.")

u = st.session_state.get("universe_ref")
if u is None:
    st.warning("Upload a trajectory and create a Universe in 📦 Trajectory Upload first.")
    st.stop()

# ---------- UI Controls ----------
colA, colB, colC, colD = st.columns([2, 1, 1, 1])

with colA:
    sel = st.text_input("MDAnalysis selection", value="nucleic")
with colB:
    max_frames = st.number_input("Frames to show (subsampled)", min_value=2, max_value=200, value=50, step=1)
with colC:
    style = st.selectbox("Style", ["cartoon", "stick", "sphere"], index=0)
with colD:
    fps = st.slider("FPS", min_value=1, max_value=30, value=8)

# Optional color mapping for common ions (hex + sphere radius)
ION_STYLE = {
    "NA+": {"hex": "#CBAACB", "radius": 0.6},
    "CL-": {"hex": "#FFFF00", "radius": 0.6},
    "K+":  {"hex": "#000000", "radius": 0.6},
    "MG":  {"hex": "#000000", "radius": 0.6},
    "CA":  {"hex": "#000000", "radius": 0.6},
    "ZN":  {"hex": "#FF0000", "radius": 1.2},
}

n_total = getattr(u.trajectory, "n_frames", 1)
if n_total < 2:
    st.info("This trajectory has a single frame. Showing a static view.")
# Decide stride to get ~max_frames models
stride = max(1, int(n_total / max(2, min(max_frames, n_total))))
st.caption(f"Total frames: **{n_total}** · Using stride: **{stride}** → up to ~**{min(max_frames, (n_total + stride - 1)//stride)}** frames shown")

# ---------- Build multi-model PDB string ----------
# We write selection for each sampled frame as a PDB and wrap with MODEL/ENDMDL
ag = u.select_atoms(sel)
models_pdb = []
tmpdir = tempfile.mkdtemp(prefix="mintllm_3d_")

# write frames to tmp and wrap with MODEL/ENDMDL
frame_index = 0
for ts in u.trajectory[::stride]:
    if frame_index >= max_frames:
        break
    out_pdb = os.path.join(tmpdir, f"frame_{frame_index}.pdb")
    ag.write(out_pdb)
    # Keep only ATOM/HETATM records (safe for py3Dmol) and wrap in MODEL blocks
    with open(out_pdb, "r") as f:
        lines = []
        for ln in f:
            if ln.startswith(("ATOM", "HETATM")) or ln.startswith(("TER", "END")):
                lines.append(ln.rstrip("\n"))
    if lines:
        block = [f"MODEL {frame_index}"] + lines + ["ENDMDL"]
        models_pdb.append("\n".join(block))
    frame_index += 1

if not models_pdb:
    st.error("No atoms in selection or failed to generate PDB frames. Try changing your selection.")
    st.stop()

multi_model_pdb = "\n".join(models_pdb)

# ---------- Build py3Dmol viewer ----------
w, h = 900, 600
view = py3Dmol.view(width=w, height=h)
# add as frames
view.addModelsAsFrames(multi_model_pdb)

# Set main style
if style == "cartoon":
    view.setStyle({'model': -1}, {'cartoon': {'color': 'spectrum'}})
elif style == "stick":
    view.setStyle({'model': -1}, {'stick': {}})
else:
    view.setStyle({'model': -1}, {'sphere': {}})

# Ion highlighting (if present)
for resn, props in ION_STYLE.items():
    view.setStyle({'model': -1, 'resn': resn}, {'sphere': {'color': props['hex'], 'radius': props['radius']}})

view.zoomTo()
# Animate — use 'interval' in ms derived from FPS
interval_ms = int(1000 / max(1, fps))
view.animate({'loop': 'forward', 'reps': 0, 'interval': interval_ms})

# ---------- Embed into Streamlit ----------
st.components.v1.html(view._make_html(), height=h, scrolling=False)

st.caption("Tip: Use 🔁 Trajectory Conversion to downsample long runs, then preview here.")
