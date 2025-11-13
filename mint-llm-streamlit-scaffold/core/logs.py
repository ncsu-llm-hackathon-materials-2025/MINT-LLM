import re
import pandas as pd

def parse_lammps_log(text: str) -> pd.DataFrame:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    headers = None
    data = []
    for ln in lines:
        if re.match(r'^\s*Step\s+', ln):
            headers = ln.split()
            data = []
            continue
        if headers and re.match(r'^\s*\d+', ln):
            row = ln.split()
            if len(row) == len(headers):
                data.append(row)
    if headers and data:
        df = pd.DataFrame(data, columns=headers).apply(pd.to_numeric, errors='ignore')
        return df
    return pd.DataFrame()

def parse_gromacs_log(text: str) -> pd.DataFrame:
    return parse_lammps_log(text)

def parse_amber_mdout(text: str) -> pd.DataFrame:
    cols = ["NSTEP","TIME(PS)","TEMP(K)","PRESS","ETOT","EPOT","EKTot"]
    rows = []
    current = {}
    for ln in text.splitlines():
        m = re.findall(r"(NSTEP|TIME\(PS\)|TEMP\(K\)|PRESS|ETOT|EPOT|EKTot)\s*=\s*([\-\d\.E\+]+)", ln)
        for k,v in m:
            current[k] = v
        if "NSTEP" in current and "TIME(PS)" in current and "TEMP(K)" in current:
            rows.append(current.copy())
            current = {}
    if rows:
        df = pd.DataFrame(rows).apply(pd.to_numeric, errors='ignore')
        return df
    return pd.DataFrame()
