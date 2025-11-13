import streamlit as st
from utils.state import init_state
from utils.ui import setup_sidebar
from core.analysis_core import compute_rdf, compute_rmsd
from core.registry import add_plot
from llm.helper import plan_from_text, chat_about_results

st.set_page_config(page_title="MINT LLM", page_icon=":herb:", layout="wide")
init_state()
setup_sidebar()

st.title("💬 Analysis Query")
st.caption("Chat with MINT about your data. Ask questions or request analyses (RDF, RMSD for now).")

# --- Guardrails: need a Universe + API key ---
api_key = st.session_state.get("openai_key", "")
if not api_key:
    st.warning("Paste your OpenAI API key in the left sidebar to chat.")
    st.stop()

u = st.session_state.get("universe_ref")
if u is None:
    st.info("Please create a Universe in 📦 Trajectory Upload first.")
    st.stop()

# --- Session chat state ---
if "chat_messages" not in st.session_state:
    st.session_state.chat_messages = []  # list of {"role": "user"/"assistant", "content": str}

# --- Optional context we pass to the LLM (small + safe) ---
context = {
    "engine": st.session_state.get("engine", ""),
    "n_atoms": getattr(getattr(u, "atoms", None), "n_atoms", None),
    "n_frames": getattr(getattr(u, "trajectory", None), "n_frames", None),
    "recent_results": [
        # You can add more compact summaries if you store them
        # For now we just reference counts; images live in the Results tab.
    ],
}

# --- Render previous messages ---
for m in st.session_state.chat_messages:
    with st.chat_message(m["role"]):
        st.markdown(m["content"])

# --- Chat input ---
prompt = st.chat_input("Ask about your data or request an analysis (e.g., RDF of water, RMSD of backbone)…")
if prompt:
    # 1) Show the user message
    st.session_state.chat_messages.append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.markdown(prompt)

    # 2) FIRST: Have the assistant respond conversationally
    with st.chat_message("assistant"):
        assistant_text = chat_about_results(prompt, api_key, context=context)
        st.markdown(assistant_text)
    st.session_state.chat_messages.append({"role": "assistant", "content": assistant_text})

    # 3) THEN: Try planning & executing analyses if the user asked for one
    # (We always try a plan; if nothing actionable is found, plan_from_text falls back safely.)
    plan = plan_from_text(prompt, api_key)
    tasks = plan.model_dump().get("tasks", [])

    ran_any = False
    for task in tasks:
        obs = task.get("observable")
        params = task.get("params", {}) or {}

        if obs == "rdf":
            defaults = {"sel_a": "all", "sel_b": "all", "rmax": 10.0, "nbins": 300}
            res = compute_rdf(u, **{**defaults, **params})
            title = f"RDF {params.get('sel_a', 'all')} vs {params.get('sel_b','all')}"
            add_plot(st.session_state, "rdf", title, params, res["png"])
            with st.chat_message("assistant"):
                st.markdown(f"✅ Ran **{title}** (peak ≈ {res['summary']['peak_r']}) and added the plot to **📈 Results**.")
            st.session_state.chat_messages.append(
                {"role": "assistant", "content": f"Ran {title} and added the plot to Results."}
            )
            ran_any = True

        elif obs == "rmsd":
            defaults = {"sel": "all", "ref_frame": 0}
            res = compute_rmsd(u, **{**defaults, **params})
            title = f"RMSD {params.get('sel','all')}"
            add_plot(st.session_state, "rmsd", title, params, res["png"])
            with st.chat_message("assistant"):
                st.markdown(
                    f"✅ Ran **{title}** (mean ≈ {res['summary']['mean']:.2f} Å) and added the plot to **📈 Results**."
                )
            st.session_state.chat_messages.append(
                {"role": "assistant", "content": f"Ran {title} and added the plot to Results."}
            )
            ran_any = True

        # (Future: elif obs in {"volume", "density"}: ...)

    if not ran_any:
        with st.chat_message("assistant"):
            st.markdown("ℹ️ No runnable analyses detected in your request. Try asking for RDF or RMSD explicitly.")
        st.session_state.chat_messages.append(
            {"role": "assistant", "content": "No runnable analyses detected. Try asking for RDF or RMSD."}
        )

    st.rerun()
