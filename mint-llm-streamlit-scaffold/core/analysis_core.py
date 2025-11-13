import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis import rms
from .plots import fig_to_png

def compute_rdf(u, sel_a='all', sel_b='all', rmax=10.0, nbins=300):
    A, B = u.select_atoms(sel_a), u.select_atoms(sel_b)
    rdf = InterRDF(A, B, range=(0.0, float(rmax)), nbins=int(nbins))
    rdf.run()
    r, g = rdf.results.bins, rdf.results.rdf
    fig, ax = plt.subplots()
    ax.plot(r, g, label='g(r)')
    ax.set_xlabel('r (Å)')
    ax.set_ylabel('g(r)')
    ax.set_title(f'RDF: {sel_a} vs {sel_b}')
    ax.legend()
    png = fig_to_png(fig)
    peak_i = int(np.argmax(g)) if len(g) else None
    summary = {'peak_r': float(r[peak_i]) if peak_i is not None else None, 'peak_g': float(g[peak_i]) if peak_i is not None else None}
    return {'r': r, 'g': g, 'png': png, 'summary': summary}

def compute_rmsd(u, sel='all', ref_frame=0, align=True):
    grp = u.select_atoms(sel)
    R = rms.RMSD(grp, grp, select=sel, ref_frame=ref_frame)
    R.run()
    arr = R.results.rmsd  # frame, time(ps or idx), rmsd(Å)
    t, y = arr[:,1], arr[:,2]
    fig, ax = plt.subplots()
    ax.plot(t, y, label='RMSD')
    ax.set_xlabel('Time (ps or frame)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title(f'RMSD: {sel}')
    ax.legend()
    png = fig_to_png(fig)
    summary = {'min': float(np.min(y)), 'max': float(np.max(y)), 'mean': float(np.mean(y))}
    return {'t': t, 'y': y, 'png': png, 'summary': summary}
