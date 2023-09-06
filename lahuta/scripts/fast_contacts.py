import unittest
import time
import tempfile

from lahuta.api import download_structures

from lahuta import Luni
from lahuta.contacts import F

RESULT = {
    "cov": (6, 2),
    "met": (0, 2),
    "carb": (30, 2),
    "hb": (0, 2),
    "whb": (0, 2),
    "ionic": (58, 2),
    "aromatic": (91, 2),
    "hydrophobic": (1203, 2),
    "phb": (863, 2),
    "wphb": (562, 2),
    "vdw": (282, 2),
    "cp": (71, 2),
    "cp2": (6, 2),
    "dp": (14, 2),
    "sp": (14, 2),
    "pp": (50, 2),
}
PDB_ID = "2RH1"


class ContactTest(unittest.TestCase):
    def test_neighbors(self) -> None:
        output_dir = tempfile.gettempdir()
        download_data = download_structures([PDB_ID], dir_loc=output_dir)

        start = time.time()
        u = Luni(download_data[PDB_ID])
        n = u.compute_neighbors(res_dif=2)

        # Compute contacAMIDE_SMARTSts
        cov = F.covalent_neighbors(n)
        self.assertEqual(cov.pairs.shape, RESULT["cov"])

        met = F.metalic_neighbors(n)
        self.assertEqual(met.pairs.shape, RESULT["met"])

        carb = F.carbonyl_neighbors(n)
        self.assertEqual(carb.pairs.shape, RESULT["carb"])

        hb = F.hbond_neighbors(n)
        self.assertEqual(hb.pairs.shape, RESULT["hb"])

        whb = F.weak_hbond_neighbors(n)
        self.assertEqual(whb.pairs.shape, RESULT["whb"])

        ionic = F.ionic_neighbors(n)
        self.assertEqual(ionic.pairs.shape, RESULT["ionic"])

        aromatic = F.aromatic_neighbors(n)
        self.assertEqual(aromatic.pairs.shape, RESULT["aromatic"])

        hydrophobic = F.hydrophobic_neighbors(n)
        self.assertEqual(hydrophobic.pairs.shape, RESULT["hydrophobic"])

        phb = F.polar_hbond_neighbors(n)
        self.assertEqual(phb.pairs.shape, RESULT["phb"])

        wphb = F.weak_polar_hbond_neighbors(n)
        self.assertEqual(wphb.pairs.shape, RESULT["wphb"])

        vdw = F.vdw_neighbors(n)
        self.assertEqual(vdw.pairs.shape, RESULT["vdw"])

        cp = F.carbon_pi(n)
        self.assertEqual(cp.pairs.shape, RESULT["cp"])

        cp2 = F.cation_pi(n)
        self.assertEqual(cp2.pairs.shape, RESULT["cp2"])

        dp = F.donor_pi(n)
        self.assertEqual(dp.pairs.shape, RESULT["dp"])

        sp = F.sulphur_pi(n)
        self.assertEqual(sp.pairs.shape, RESULT["sp"])

        pp = F.plane_plane_neighbors(n)
        self.assertEqual(pp.pairs.shape, RESULT["pp"])

        end = time.time()

        self.assertEqual(aromatic.pairs.shape, (91, 2))
        print(f"Time elapsed: {end - start}")


if __name__ == "__main__":
    unittest.main()
