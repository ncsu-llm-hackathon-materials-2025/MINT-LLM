import streamlit as st
import base64, io, zipfile, json
from utils.state import init_state
from utils.ui import setup_sidebar

st.set_page_config(page_title="MINT LLM", page_icon=":herb:", layout="wide")
init_state()
setup_sidebar()

st.title("📈 Results")
st.caption("All generated plots and tables appear here. Download individually or as a bundle.")

items = st.session_state.get("results_registry", [])
if not items:
    st.info("No results yet. Run analyses in Assessment or Analysis Query.")
else:
    cols = st.columns(3)
    for i, it in enumerate(items):
        with cols[i % 3]:
            st.subheader(it["title"])
            st.image(base64.b64decode(it["png_b64"]), use_column_width=True)
            st.caption(f"Kind: {it['kind']} • Created: {it['created_at']}")
            if it.get("csv_b64"):
                st.download_button("Download CSV", base64.b64decode(it["csv_b64"]), file_name=f"{it['id']}.csv", mime="text/csv")

    if st.button("Download All (zip)"):
        mem = io.BytesIO()
        with zipfile.ZipFile(mem, mode="w") as z:
            for it in items:
                z.writestr(f"{it['id']}_{it['title']}.png", base64.b64decode(it["png_b64"]))
                if it.get("csv_b64"):
                    z.writestr(f"{it['id']}_{it['title']}.csv", base64.b64decode(it["csv_b64"]))
            z.writestr("summary.json", json.dumps(items, indent=2))
        st.download_button("Save Zip", mem.getvalue(), file_name="mint-llm-results.zip", mime="application/zip")
