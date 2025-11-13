import streamlit as st

DEFAULTS = {
    "openai_key": "",
    "engine": "",
    "uploads": {},              # { 'lammps': {...}, 'gromacs': {...}, 'amber': {...} }
    "universe_ref": None,       # placeholder to cache Universe handle
    "results_registry": [],     # list of artifacts
    "report_blocks": [],
}

def init_state():
    for k, v in DEFAULTS.items():
        if k not in st.session_state:
            st.session_state[k] = v
