import time, base64

def add_plot(ss, kind: str, title: str, params: dict, png_bytes: bytes, csv_bytes: bytes|None=None, tags: list[str]|None=None):
    item = {
        "id": f"{int(time.time()*1000)}",
        "kind": kind,
        "title": title,
        "params": params,
        "png_b64": base64.b64encode(png_bytes).decode(),
        "csv_b64": base64.b64encode(csv_bytes).decode() if csv_bytes else None,
        "tags": tags or [],
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    ss["results_registry"].append(item)

def list_items(ss, kind: str|None=None):
    items = ss.get("results_registry", [])
    if kind:
        items = [i for i in items if i["kind"] == kind]
    return list(reversed(items))
