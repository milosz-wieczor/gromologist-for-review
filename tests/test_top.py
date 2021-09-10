import gromologist as gml
import unittest


class BasicTopTest(unittest.TestCase):

    def setUp(self) -> None:
        self.top = gml.Top('pentapeptide.top')

    def tearDown(self) -> None:
        del self.top

    def test_natoms(self):
        nat = len(self.top.atoms)
        self.assertEqual(nat, 87)

    def test_atom_charge(self):
        charge = self.top.atoms[0].charge
        self.assertEqual(charge, -0.3)

    def test_nmols(self):
        nmol = len(self.top.molecules)
        self.assertEqual(nmol, 10)

    def test_clear_mols(self):
        self.top.clear_sections()
        nmol = len(self.top.molecules)
        self.assertEqual(nmol, 1)

    def test_del_atom(self):
        pentapeptide = self.top.molecules[0]
        pentapeptide.del_atom(1)
        nat = len(self.top.atoms)
        self.assertEqual(nat, 86)

    def test_add_params_bond(self):
        self.top.add_ff_params()
        self.assertEqual(len(self.top.molecules[0].subsections[2].entries[3].params_state_a), 2)

    def test_add_params_angle(self):
        self.top.add_ff_params()
        self.assertEqual(len(self.top.molecules[0].subsections[4].entries[3].params_state_a), 4)


if __name__ == "__main__":
    unittest.main()