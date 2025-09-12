# file: io_convert.py
from pathlib import Path
import sys
import numpy as np
import MDAnalysis as mda

# !-------------------------------------------------------------------------------------- !
READABLE = {".pdb", ".xyz", ".xtc", ".trr", ".dcd",".lammpstrj"}
WRITABLE = {".pdb", ".xyz", ".xtc", ".trr", ".dcd", ".gro", ".lammpstrj"}
# !-------------------------------------------------------------------------------------- !

#===================================== Helper Functions ============================================
# trim suffix from the filename and return file extension in lowercase (eg. ".pdb")
def _ext(p):
    return Path(p).suffix.lower()

# requires pdb for coordinate-only formats
def NeedTopo(infile):
    return _ext(infile) in {".xtc", ".trr", ".dcd", ".lammpstrj", ".xyz"}

# validate periodicity [A,B,C,alpha,beta,gamma]
def ValidDims(d):
    try:
        d = np.asarray(d, float)
        return d.shape[0] == 6 and np.all(np.isfinite(d)) and np.all(d[:3] > 0.0)
    except Exception:
        return False

# pull box from topology file: try MDAnalysis else fall back to ASE for PDB reader
def DimsFromTopology(topology):
    if not topology:
        return None
    # Try MDAnalysis : (works for PDB/GRO/TPR that carry a box)
    try:
        ut = mda.Universe(topology)
        d = getattr(ut.trajectory.ts, "dimensions", None)
        if ValidDims(d):
            return np.array(d, float)
    except Exception:
        pass
    # Fall back to ASE PDB reader
    if _ext(topology) == ".pdb":
        try:
            from ase.io import read as ase_read
            at = ase_read(topology)
            L = np.array(at.get_cell().lengths(), float)
            A = np.array(at.get_cell().angles(),  float)
            if np.all(L > 0.0):
                return np.array([L[0], L[1], L[2], A[0], A[1], A[2]], float)
        except Exception:
            pass
    return None

#=====================================================================================================
# standard MDAnalysis function to read/write various formats
def MdaUniverse(infile, topology):
    """
        infile: str, topology: Optional[str]
        Standard MDAnalysis Universe open
    """
    if NeedTopo(infile) and topology is None:
        raise ValueError(f"{_ext(infile)} is coordinate-only. Please provide a topology (.tpr/.gro/.pdb/.prmtop).")
    return mda.Universe(topology, infile) if topology else mda.Universe(infile)

#========================================================================================================
# Speacial function to read LAMMPS trajectory with ASE
# and restore atomic numbers + cell from a PDB file.
# This is needed because MDAnalysis cannot read atomic numbers from LAMMPS dump files,
# and LAMMPS dump files do not contain element symbols.
def PDBWrapper(lmp_file, pdb_file):
    """
    lmp_file: str, pdb_file: str
    Read '.lammpstrj' using ASE, copy atomic numbers & cell from PDB file, then build a MDAnalysis Universe for all frames.
    ```python
        >>> import MDAnalysis as mda
        >>> u = PDBWrapper("test.lammpstrj", "peptides.pdb")
        >>> with mda.Writer("out.xtc", n_atoms=u.atoms.n_atoms, multiframe=True) as w:
        >>>     for ts in u.trajectory:
        >>>         w.write(u.atoms)
    ```
    """
    from ase.io import read as ase_read

    # Topology of atoms from PDB file using ASE (atomic numbers, lattice vectors)
    top_atoms = ase_read(pdb_file)
    Z = np.asarray(top_atoms.get_atomic_numbers(), int)     # get atomic numbers
    cell = top_atoms.get_cell()                             # get lattice vectors [A,B,C,alpha,beta,gamma]
    L = np.array(cell.lengths(), float)                     # [ABC]
    A = np.array(cell.angles(), float)                      # [alpha,beta,gamma]
    if not np.all(L > 0):
        msg = "\n".join([
            f"PDB '{pdb_file}' has no valid periodicity or lattice parameters (CRYST1 missing or zero lengths).",
            f"Hint: see for line on top of the '{pdb_file}' file.",
            "CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1", ])
        raise ValueError(msg)

    # Read all frames using ASE read
    images = ase_read(lmp_file, index=":")  # list[Atoms] (even if 1 frame)
    if not isinstance(images, list):
        images = [images]
    n_frames = len(images)
    n_atoms = len(top_atoms)

    # Build empty list XYZ and per-frame dimensions [A,B,C,alpha,beta,gamma]
    xyz = np.empty((n_frames, n_atoms, 3), dtype=float)
    dims = np.tile(np.array([L[0], L[1], L[2], A[0], A[1], A[2]], float), (n_frames, 1))

    for i, a in enumerate(images):
        if len(a) != n_atoms:
            raise ValueError(f"Atom count mismatch in frame {i}: dump={len(a)}, PDB={n_atoms}")
        # Overwrite atomic numbers and cell from PDB
        a.set_atomic_numbers(Z)
        a.set_cell(cell, scale_atoms=False)
        a.set_pbc(True)
        xyz[i] = a.get_positions()

    # Build a MDAnalysis Universe for all frames
    u = mda.Universe(pdb_file)
    u.load_new(xyz, dimensions=dims)
    return u
#========================================================================================================
# returns an MDAnalysis Universe for trajectory conversion
def MdanalysisWriter(u, outfile, xtc_precision=3):
    """
        u: mda.Universe, outfile: str, xtc_precision: int = 3
    """
    e = _ext(outfile)
    if e in {".pdb", ".xyz"}:
        return mda.Writer(str(outfile), multiframe=True)
    if e == ".xtc":
        return mda.Writer(str(outfile), n_atoms=u.atoms.n_atoms, multiframe=True, precision=xtc_precision)
    if e in {".trr", ".dcd"}:
        return mda.Writer(str(outfile), n_atoms=u.atoms.n_atoms, multiframe=True)
    return None
#========================================================================================================
# Spacial: only write a single frame to GRO format
def WriteGRO(u, outfile, frame_index):
    """
        u: mda.Universe, outfile: str, frame_index: int
        GROMACS `gro` file writer: single frame only.
        Coordinates are written in nm after converting from agstrom, and orthorhombic box if present.
    ```python
        >>> from io_trajcon import TrajConvert
        >>> TrajConvert("traj.xyz", "snap.gro", topology="peptides.pdb", start=64)
        >>>     [GRO] Writing single frame -> 'snap.gro'  (frame=64, atoms=82)
    ```
    """
    ts = u.trajectory[frame_index]
    nat = u.atoms.n_atoms
    pos_nm = ts.positions / 10.0  # angstrom -> nm
    print(f"[GRO] Writing single frame -> '{outfile}'  (frame={frame_index}, atoms={nat})")
    with open(outfile, "w") as f:
        title = getattr(u, "filename", "Converted by io_trajcon.py")
        f.write(f"{title}\n")
        f.write(f"{nat:5d}\n")

        resids = u.atoms.resids
        resnames = u.atoms.resnames
        names = u.atoms.names
        for i in range(nat):
            resid = int(resids[i]) % 100000
            resn = (resnames[i] or "")[:5]
            name = (names[i] or "")[:5]
            atomnr = (i + 1) % 100000
            x, y, z = pos_nm[i]
            f.write(f"{resid:5d}{resn:<5s}{name:>5s}{atomnr:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")

        dims = getattr(ts, "dimensions", None)
        if dims is not None and len(dims) == 6 and np.isclose(dims[3:], [90, 90, 90]).all():
            A, B, C = dims[:3]
            f.write(f"{A/10.0:10.5f} {B/10.0:10.5f} {C/10.0:10.5f}\n")
        else:
            f.write("0.00000 0.00000 0.00000\n")
            print("[WARN] Writing GRO with zero box (no valid orthorhombic box found).", file=sys.stderr)
#========================================================================================================
# extract per-element integer 'type' from PDB elements using ASE for lammpstrj output
def TypesExtractASE(pdb_file):
    """
        Build per-element integer 'type' list from PDB file using ASE.
        Returns (types_list, atoms) where atoms is the ASE Atoms object (for cell info).
        - sorted by atomic symbol: ['H','C','N','O']
        - assign a integer type : {'H':1, 'C':2, 'N':3, 'O':4}
        
    ``` python
    >>> from io_trajcon import TypesExtractASE
    >>> type, atoms = TypesExtractASE("topology.pdb")
    >>> Elements found: [H, C, N, O] -> mapped to types: {H:1, C:2, N:3, O:4}
    ```
    """
    from ase.io import read as ase_read
    from ase.data import chemical_symbols as _cs
    atoms = ase_read(pdb_file)
    syms = atoms.get_chemical_symbols()
    def _key(x):
        try:
            return (_cs.index(x), x)
        except ValueError:
            return (999, x)
    uniq = sorted(set(syms), key=_key)
    lut = {el: i + 1 for i, el in enumerate(uniq)}
    types_list = [lut[s] for s in syms]
    print(f"Elements found: [{', '.join(sorted(lut, key=lut.get))}] -> mapped to types: " 
          + "{" + ", ".join(f"{el}:{lut[el]}" for el in sorted(lut, key=lut.get)) + "}")
    return types_list, atoms
#========================================================================================================
# write LAMMPS lammpstrj
def WriteLAMMPStrj(u , outfile, frames, pdb_file):
    """
        - u : mda.Universe , outfile : str, frames : str, pdb_file : str
        Peridic box from PDB CRYST1 and wrapped positions to the center of the box.
    """
    types_list, top_atoms = TypesExtractASE(pdb_file)
    if len(types_list) != u.atoms.n_atoms:
        raise ValueError("Atom count mismatch between PDB and Trajectory file.")
    L_A = np.array(top_atoms.get_cell().lengths(), float)    # ABC in angstrom
    if not np.all(L_A > 0):
        raise ValueError("PDB file has no valid lattice information.")
    top_atoms.set_pbc(True)
    nat = u.atoms.n_atoms

    # LAMMPS lammpstrj style formatting
    with open(outfile, "w") as f:
        for i in frames:
            ts = u.trajectory[i]                # atomic positions
            fr = top_atoms.copy()               # create a copy of PDB atoms
            fr.set_positions(ts.positions)      # Overlay PDB to the trajectory
            fr.set_pbc(True)                    # Set periodicity if not found
            fr.wrap()                           # wrap atom inside the box
            pos = fr.get_positions()            # Å
            pos -= pos.min(axis=0)              # fix origin to 0

            f.write("ITEM: TIMESTEP\n")
            f.write(f"{ts.frame}\n")
            f.write("ITEM: NUMBER OF ATOMS\n")
            f.write(f"{nat}\n")
            f.write("ITEM: BOX BOUNDS pp pp pp\n")
            f.write(f"0.0 {L_A[0]:.10f}\n0.0 {L_A[1]:.10f}\n0.0 {L_A[2]:.10f}\n")
            f.write("ITEM: ATOMS id type x y z\n")
            for idx in range(nat):
                x, y, z = pos[idx]
                f.write(f"{idx+1} {types_list[idx]} {x:.10f} {y:.10f} {z:.10f}\n")
#========================================================================================================
# Main entry point for format conversion
def TrajConvert(infile, outfile, *, topology = None, start = None, stop = None, step = None, xtc_precision = 3):
    """
    Convert between PDB, XYZ, XTC, TRR, DCD, GRO, lammpstrj.
    """
    in_ext, out_ext = _ext(infile), _ext(outfile)
    print(f"[IO] InTraj: '{infile}' [{in_ext}]  =>  OutTraj: '{outfile}' [{out_ext}]")
    print(f"[IO] Topology: {repr(topology) if topology else '(none)'}  |  Slice: start={start}, stop={stop}, step={step}")
    if in_ext not in READABLE:
        raise ValueError(f"Unsupported input: {in_ext}")
    if out_ext not in WRITABLE:
        if out_ext == ".tpr":
            raise ValueError("Unsupported output: .tpr cannot be written; use .gro/.pdb/.xtc/.trr instead.")
        raise ValueError(f"Unsupported output: {out_ext}")

    # Build Universe (special case for reading lammpstrj)
    if in_ext == ".lammpstrj":
        if not topology or _ext(topology) != ".pdb":
            raise ValueError("Reading .lammpstrj requires a PDB topology to restore atomic symbols.")
        print("[READ] Using PDBWrapper (LAMMPS dump + PDB for atomic_symbols/box).")
        u = PDBWrapper(infile, topology)
    else:
        if NeedTopo(infile) and not topology:
            print("[WARN] Coord-only input without topology will raise an error.", file=sys.stderr)
        u = MdaUniverse(infile, topology)

    frame_slice = slice(start, stop, step)
    frames = range(u.trajectory.n_frames)[frame_slice]

    # Fetch box once from topology file
    dims_from_top = DimsFromTopology(topology)

    # Non-LAMMPS outputs using MDAnalysis writers
    W = MdanalysisWriter(u, outfile, xtc_precision)
    if W is not None:
        with W as w:
            for i in frames:
                ts = u.trajectory[i]
                # inject box if frame lacks valid dimensions
                d = getattr(ts, "dimensions", None)
                if (not ValidDims(d)) and ValidDims(dims_from_top):
                    ts.dimensions = dims_from_top
                w.write(u.atoms)
        return

    # GRO (single-frame)
    if out_ext == ".gro":
        try:
            first_frame = next(iter(frames))
        except StopIteration:
            first_frame = 0
        #
        ts0 = u.trajectory[first_frame]
        d0 = getattr(ts0, "dimensions", None)
        if (not ValidDims(d0)) and ValidDims(dims_from_top):
            ts0.dimensions = dims_from_top
        WriteGRO(u, outfile, first_frame)
        return

    # LAMMPS lammpstrj format requires box information from a topology file (eg: '.pdb')
    if out_ext == ".lammpstrj":
        if not topology or _ext(topology) != ".pdb":
            raise ValueError("LAMMPS output requires a PDB topology and element types (eg: '.pdb').")
        WriteLAMMPStrj(u, outfile, frames, topology)
        return

    raise ValueError("Unhandled output format.")
#========================================================================================================