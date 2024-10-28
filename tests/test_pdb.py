import gromologist as gml
import unittest
import os
import numpy as np


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
        [self.pdb.delete_atom(1, renumber=True) for _ in range(6)]
        self.assertEqual(self.pdb.natoms, 80)

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

    def test_altlocs(self):
        self.assertEqual(self.pdb.altlocs, ' ')

    def test_conect(self):
        self.pdb.add_conect()
        self.assertEqual(self.pdb.conect[5], [1, 7, 11, 6])

    def test_chains(self):
        self.assertEqual(self.pdb.chains, [' '])
        self.pdb.add_chains()
        self.assertEqual(self.pdb.chains, ['A'])

    def test_addqt(self):
        self.pdb.add_top('pentapeptide.top')
        self.pdb.check_top()
        self.assertEqual(self.pdb.qt, False)
        self.pdb.add_qt()
        self.assertEqual(self.pdb.atoms[11].q, -0.51)
        self.assertEqual(self.pdb.qt, True)

    def test_transl(self):
        self.assertEqual(self.pdb.get_coords()[0,0] * self.pdb.get_coords()[20,2], -27.0584)
        self.pdb.translate_selection('resid 1', np.array([1, 0, 1]))
        self.assertEqual(self.pdb.get_coords()[0, 0] * self.pdb.get_coords()[20, 2], -17.9784)

    def test_box(self):
        self.assertEqual(np.sum(self.pdb.box), 556.828)
        self.pdb.save_gro('test.gro')
        p = gml.Pdb('test.gro')
        self.assertEqual(np.sum(p.box), 556.828)
        del p

    def test_beta(self):
        self.pdb.set_beta(np.arange(2, 7), 'name CA')
        self.assertEqual(self.pdb.atoms[14].beta, 3)
        self.pdb.set_beta(np.arange(2, 7), 'name CA', smooth=15)
        self.assertEqual(self.pdb.atoms[14].beta, 3.192584474529868)

    def test_elements_coords(self):
        self.assertEqual(self.pdb.get_coords().shape, (87, 3))
        self.assertEqual(self.pdb.get_coords('element H').shape, (0,))
        self.pdb.add_elements()
        self.assertEqual(self.pdb.get_coords('element H').shape, (47, 3))

    def test_residue(self):
        self.assertEqual(len(self.pdb.residues), 5)
        self.assertEqual(self.pdb.residues[2].selection, 'resid 3 and resname PHE and chain A')

    def test_protein_seq(self):
        self.assertEqual(self.pdb.print_protein_sequence(), [''])
        self.pdb.add_chains()
        self.assertEqual(self.pdb.print_protein_sequence()[0], 'ALFIV')


if __name__ == "__main__":
    unittest.main()

