import streamlit as st
from utils.state import init_state
from utils.ui import setup_sidebar

st.set_page_config(page_title="MINT LLM", page_icon=":herb:", layout="wide")
init_state()
setup_sidebar()

# st.set_page_config(
#     page_title="MINT LLM",
#     page_icon=":herb:",
#     layout="wide",
# )

init_state()

st.title("MINT LLM — Molecular INsights Toolkit")
st.caption("Unified, LLM-assisted analysis for soft-matter molecular dynamics")

# # Sidebar: global controls
# with st.sidebar:
#     st.header("Global Settings")
#     st.text_input("OpenAI API Key", key="openai_key", type="password", help="Used only for this session.")
#     st.selectbox("Select Engine", ["", "LAMMPS", "GROMACS", "AMBER"], key="engine")
#     st.divider()
#     st.markdown("**Quick Links**")
#     st.page_link("pages/1_Simulation_Assessment.py", label="🧪 Simulation Assessment")
#     st.page_link("pages/2_Trajectory_Upload.py", label="📦 Trajectory Upload")
#     st.page_link("pages/3_Analysis_Query.py", label="💬 Analysis Query")
#     st.page_link("pages/4_Results.py", label="📈 Results")
#     st.page_link("pages/5_3D_Viewer.py", label="🧬 3D Viewer")
#     st.page_link("pages/6_Trajectory_Conversion.py", label="🔁 Trajectory Conversion")
#     if st.button("Reset Session", type="primary"):
#         for k in list(st.session_state.keys()):
#             del st.session_state[k]
#         st.rerun()

st.markdown(
    '''
    ### Welcome
    MINT LLM helps you assess simulations from logs, upload trajectories for analysis (RDF, RMSD, volume, density), chat with an LLM about your data, and collect results into a report.

    Workflow: Home -> Assessment -> Upload -> Query -> Results -> 3D Viewer -> Conversion
    '''
)
st.info("Tip: Set your engine and paste your OpenAI API key in the sidebar first.")
