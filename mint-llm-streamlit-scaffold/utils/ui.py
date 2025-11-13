# utils/ui.py
import streamlit as st

ENGINE_OPTIONS = ["", "LAMMPS", "GROMACS", "AMBER"]

def setup_sidebar():
    with st.sidebar:
        st.header("Global Settings")

        # --- Canonical state defaults (persist across pages) ---
        if "openai_key" not in st.session_state:
            st.session_state.openai_key = ""
        if "engine" not in st.session_state:
            st.session_state.engine = ""

        # --- Engine: use a proxy widget key and sync back to canonical ---
        eng_index = ENGINE_OPTIONS.index(st.session_state.engine) if st.session_state.engine in ENGINE_OPTIONS else 0
        engine_widget = st.selectbox(
            "Select Engine",
            ENGINE_OPTIONS,
            index=eng_index,
            key="__engine_select",   # proxy key
            help="This choice is shared across all pages."
        )
        # Sync to canonical state
        st.session_state.engine = engine_widget

        # --- API key: password input is intentionally blank on rerun; we keep canonical in session ---
        entered_key = st.text_input(
            "OpenAI API Key",
            type="password",
            key="__openai_key_input",  # proxy key
            help="Stored only in memory for this session."
        )
        if entered_key:
            st.session_state.openai_key = entered_key

        # Status indicators so you can *see* what's set even if the field looks empty
        if st.session_state.openai_key:
            st.success("API key is set for this session")
        else:
            st.warning("API key not set")

        if st.session_state.engine:
            st.caption(f"Engine: **{st.session_state.engine}**")
        else:
            st.caption("Engine: *(not selected)*")

        st.divider()
        if st.button("Reset Session", type="primary"):
            # Clear canonical & any page-local keys; then rerun
            for k in list(st.session_state.keys()):
                del st.session_state[k]
            st.rerun()
