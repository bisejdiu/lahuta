import unittest
import tempfile
import time
from lahuta.api import download_structures
from lahuta import Luni
from lahuta.contacts import F


class SimpleTest(unittest.TestCase):
    def test_aromatic_neighbors(self) -> None:
        output_dir = tempfile.gettempdir()
        pdb_id = "2RH1"
        download_data = download_structures([pdb_id], dir_loc=output_dir)

        start = time.time()
        u = Luni(download_data[pdb_id])
        n = u.compute_neighbors(res_dif=2)
        aromatic = F.aromatic_neighbors(n)
        end = time.time()

        self.assertEqual(aromatic.pairs.shape, (91, 2))
        print(f"Time elapsed: {end - start}")


if __name__ == "__main__":
    unittest.main()
