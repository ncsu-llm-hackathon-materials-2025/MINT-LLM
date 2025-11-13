import os, uuid, pathlib

def save_uploaded_file(uploaded_file, subdir="uploads") -> str:
    base = pathlib.Path(".cache")
    base.mkdir(exist_ok=True)
    (base / subdir).mkdir(parents=True, exist_ok=True)
    suffix = os.path.splitext(uploaded_file.name)[1]
    filename = f"{uuid.uuid4().hex}{suffix}"
    path = base / subdir / filename
    with open(path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return str(path)

def validate_ext(name: str, allowed: tuple[str, ...]) -> bool:
    return name.lower().endswith(allowed)
