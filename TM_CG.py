import mdtraj as md
import numpy as np
import pandas as pd

def cg_frame2frame_rmsd(topology_file, trajectory_file, stride=1, pdb_out="aligned.pdb", csv_out="frame_to_frame_rmsd.csv"):
    """
    Compute per-residue frame-to-frame RMSD for a coarse-grained trajectory.

    Parameters
    ----------
    topology_file : str
        Path to the topology file (e.g., .gro)
    trajectory_file : str
        Path to the trajectory file (e.g., .xtc)
    stride : int
        Load every nth frame
    pdb_out : str
        File name to save aligned PDB
    csv_out : str
        File name to save RMSD CSV
    """

    print("Loading trajectory...")
    traj = md.load(trajectory_file, top=topology_file, stride=stride)
    n_residues = traj.n_residues
    n_frames = traj.n_frames
    print(f"Residues: {n_residues}, Frames: {n_frames}")

    # Align all frames to first frame
    traj.superpose(traj, 0)
    traj.save_pdb(pdb_out)
    print(f"Aligned trajectory saved as: {pdb_out}")

    # Initialize RMSD matrix: residues x (frames-1)
    rmsd_matrix = np.zeros((n_residues, n_frames-1))

    print("Computing frame-to-frame RMSD per residue...")
    for i, res in enumerate(traj.topology.residues):
        atom_indices = [atom.index for atom in res.atoms]  # usually 1 backbone bead per residue in CG
        for t in range(1, n_frames):
            diff = traj.xyz[t, atom_indices, :] - traj.xyz[t-1, atom_indices, :]
            rmsd_nm = np.sqrt(np.mean(np.sum(diff**2, axis=1)))  # RMSD in nm
            rmsd_matrix[i, t-1] = rmsd_nm

    # Optional: convert to Å
    rmsd_matrix_angstrom = rmsd_matrix * 10

    # Save CSV
    pd.DataFrame(rmsd_matrix_angstrom).to_csv(csv_out, index=False)
    print(f"Frame-to-frame RMSD saved as: {csv_out}")

    return rmsd_matrix_angstrom

# ======================
# Example usage
# ======================
if __name__ == "__main__":
    topology_file = "12_md20_2_0.67_SpoIIIEinserted_2A_10ns_R2centered.gro"       # replace with your PDB/topology
    trajectory_file = "12_md20_2_0.67_SpoIIIEinserted_2A_10ns_R2centered.xtc" # replace with your trajectory
    rmsd_matrix = cg_frame2frame_rmsd(
        topology_file,
        trajectory_file,
        stride=1,
        pdb_out="aligned.pdb",
        csv_out="12_SpoIIIEmemins_RMSD.csv"
    )
