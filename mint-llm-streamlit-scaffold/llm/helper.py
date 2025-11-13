from typing import Dict, Any, List
from pydantic import BaseModel, Field, ValidationError
import json
import openai

class Task(BaseModel):
    observable: str
    params: Dict[str, Any] = Field(default_factory=dict)

class Plan(BaseModel):
    tasks: List[Task] = Field(default_factory=list)

SYSTEM_PROMPT = (
    "You are MINT LLM, an MD analysis planner. "
    "Return a strict JSON object with 'tasks', each having 'observable' in ['rdf','rmsd','volume','density'] and 'params'. "
    "Default params: rdf: {sel_a:'all', sel_b:'all', rmax:10.0, nbins:300}; rmsd: {sel:'all', ref_frame:0}; "
    "Return ONLY JSON."
)

def get_client(api_key: str):
    client = openai.OpenAI(api_key=api_key)
    return client

def plan_from_text(prompt: str, api_key: str) -> Plan:
    client = get_client(api_key)
    resp = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {"role":"system","content":SYSTEM_PROMPT},
            {"role":"user","content":prompt}
        ],
        temperature=0.2,
    )
    content = resp.choices[0].message.content
    try:
        data = json.loads(content)
        return Plan(**data)
    except (json.JSONDecodeError, ValidationError):
        return Plan(tasks=[Task(observable="rdf", params={"sel_a":"all","sel_b":"all","rmax":10.0,"nbins":300})])

# --- Simple chat about results / data ---
def chat_about_results(user_text: str, api_key: str, context: dict | None = None) -> str:
    """
    Send a normal chat message to the model with an optional context dict
    describing the current dataset/results (engine, selections, stats).
    Returns assistant text.
    """
    client = get_client(api_key)
    sys = (
        "You are MINT LLM, a helpful molecular dynamics analysis assistant. "
        "Be concise, cautious with assumptions, and suggest follow-up analyses when useful. "
        "If asked to run an analysis, describe it briefly; a separate tool will execute it."
    )
    ctx_txt = ""
    if context:
        import json
        ctx_txt = "CONTEXT (JSON):\n" + json.dumps(context, indent=2)

    resp = client.chat.completions.create(
        model="gpt-4o-mini",
        temperature=0.3,
        messages=[
            {"role": "system", "content": sys},
            {"role": "user", "content": f"{user_text}\n\n{ctx_txt}".strip()},
        ],
    )
    return resp.choices[0].message.content

