import MDAnalysis as mda
from MDAnalysis.analysis import rms

def run_rmsd(topology_file, trajectory_file):
    u = mda.Universe(topology_file, trajectory_file)
    R = rms.RMSD(u, u, select="backbone").run()
    return R.rmsd.tolist()

