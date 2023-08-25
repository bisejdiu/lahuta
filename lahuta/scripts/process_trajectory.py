from pathlib import Path

import MDAnalysis as mda

from lahuta import Luni
from lahuta.contacts.computer import LahutaContacts, LahutaTrajectoryContacts


def remove_hidden_files(test_dir: Path) -> None:
    """Remove hidden files that are created by the tests.

    Args:
        test_dir (Path): Path to the test directory.
    """

    # Find hidden files that end with .lock or .npz in that directory
    for pattern in ('.*.lock', '.*.npz'):
        hidden_files = test_dir.glob(pattern)
        for file_path in hidden_files:
            if file_path.is_file():
                file_path.unlink()


if __name__ == "__main__":
    file_dir = Path(__file__).parent.parent / "tests" / "data"
    # Load the universe
    coords = str(file_dir / 'conf0_197_sel.pdb')
    traj = str(file_dir / 'trj_197_sel.xtc')

    selection = 'protein'
    u = mda.Universe(coords, traj).select_atoms(selection)

    luni = Luni(u.atoms)
    luni.ready()

    ct = LahutaContacts(contact_type='atom-atom')
    trj_contacts = LahutaTrajectoryContacts(res_dif=2, radius=5)

    trj_contacts.compute(luni, ct, n_jobs=4)

    print(trj_contacts.results.get(0))

    remove_hidden_files(file_dir)
