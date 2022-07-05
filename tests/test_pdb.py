import gromologist as gml
import unittest
import os


class BasicTopTest(unittest.TestCase):

    def setUp(self) -> None:
        self.pdb = gml.Pdb('pentapeptide.pdb')

    def tearDown(self) -> None:
        del self.pdb

    def test_natoms(self):
        # checks the pdb.natoms attr
        self.assertEqual(self.pdb.natoms, 87)

    def test_del_atom(self):
        # checks if deleting atoms works
        self.pdb.delete_atom(11)
        self.assertEqual(self.pdb.natoms, 86)

    def test_add_atom(self):
        # checks if adding atoms works
        self.pdb.insert_atom(20, 'DUM', hooksel='name O and resid 2', bondlength=1.1, p1_sel='name CA and resid 1',
                             p2_sel='name CA and resid 2')
        self.assertEqual(self.pdb.natoms, 88)

    def test_ala_mut(self):
        # checks if atoms are removed/added correctly by the mutation module
        self.pdb.mutate_protein_residue(3, 'A')
        self.assertEqual(self.pdb.natoms, 77)

    def test_selection(self):
        # checks if complex selection expressions work
        self.assertEqual(len(self.pdb.get_atoms('same residue as within 3 of serial 20')), 19)


if __name__ == "__main__":
    unittest.main()

