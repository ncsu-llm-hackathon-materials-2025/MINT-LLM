import MDAnalysis as mda
from MDAnalysis.analysis import rms

def run_rmsd(topology_file, trajectory_file):
    print("using the provided RMSD function")
    u = mda.Universe(topology_file, trajectory_file)
    R = rms.RMSD(u, u, select="backbone").run()
    return R.rmsd.tolist()

