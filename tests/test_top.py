import gromologist as gml
import unittest


class BasicTopTest(unittest.TestCase):

    def setUp(self) -> None:
        self.top = gml.Top('pentapeptide.top')

    def tearDown(self) -> None:
        del self.top

    def test_natoms(self):
        # checks the top.atoms list
        nat = len(self.top.atoms)
        self.assertEqual(nat, 87)

    def test_atom_charge(self):
        # checks if charges are summed correctly
        charge = self.top.atoms[0].charge
        self.assertEqual(charge, -0.3)

    def test_nmols(self):
        # checks if molecules are counted correctly
        nmol = len(self.top.molecules)
        self.assertEqual(nmol, 10)

    def test_clear_mols(self):
        # checks cleaning up unused molecule definitions
        self.top.clear_sections()
        nmol = len(self.top.molecules)
        self.assertEqual(nmol, 1)

    def test_del_atom(self):
        # checks if deleting atoms works
        pentapeptide = self.top.molecules[0]
        pentapeptide.del_atom(1)
        nat = len(self.top.atoms)
        self.assertEqual(nat, 86)

    def test_add_params_bond(self):
        # checks if explicit adding of FF parameters for bonds works
        self.top.add_ff_params()
        self.assertEqual(len(self.top.molecules[0].subsections[2].entries[3].params_state_a), 2)

    def test_add_params_angle(self):
        # checks if explicit adding of FF parameters for angles works
        self.top.add_ff_params()
        self.assertEqual(len(self.top.molecules[0].subsections[4].entries[3].params_state_a), 4)


if __name__ == "__main__":
    unittest.main()