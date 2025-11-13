import MDAnalysis as mda

# def load_universe(engine: str, files: dict) -> mda.Universe:
#     e = (engine or '').lower()
#     if e == 'lammps':
#         return mda.Universe(files['data.lammps'], files['dump.lammpstrj'][0])
#     if e == 'gromacs':
#         return mda.Universe(files['topol.tpr'], files['traj'][0])
#     if e == 'amber':
#         return mda.Universe(files['prmtop'], files['prod'][0])
#     raise ValueError('Unsupported engine or missing files')

def load_universe(engine: str, files: dict) -> mda.Universe:
    e = (engine or '').lower()
    if e == 'lammps':
        return mda.Universe(files['data.lammps'], files['dump.lammpstrj'][0])
    if e == 'gromacs':
        # Accept either TPR or GRO as topology
        if 'topol.tpr' in files:
            topo = files['topol.tpr']
        elif 'topol.gro' in files:
            topo = files['topol.gro']
        else:
            raise ValueError('GROMACS topology not provided (need .tpr or .gro)')
        return mda.Universe(topo, files['traj'][0])
    if e == 'amber':
        return mda.Universe(files['prmtop'], files['prod'][0])
    raise ValueError('Unsupported engine or missing files')
