import gromologist as gml
import unittest
import os


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

    def test_add_atom(self):
        # checks if adding atoms works
        pentapeptide = self.top.molecules[0]
        pentapeptide.add_atom(20, 'DUM', 'NA', charge=-1.0, resid=4, resname='CHK', mass=2.5)
        nat = len(self.top.atoms)
        self.assertEqual(nat, 88)

    def test_add_atom_charge(self):
        # checks if adding atoms updates total charge
        pentapeptide = self.top.molecules[0]
        orig_charge = self.top.charge
        pentapeptide.add_atom(20, 'DUM', 'NA', charge=-1.0, resid=4, resname='CHK', mass=2.5)
        self.assertAlmostEqual(self.top.charge, orig_charge-1)

    def test_ala_mut(self):
        # checks if atoms are removed/added correctly by the mutation module
        pentapeptide = self.top.molecules[0]
        pentapeptide.mutate_protein_residue(3, 'A', rtp='merged.rtp')
        self.assertEqual(self.top.natoms, 77)

    def test_add_params_bond(self):
        # checks if explicit adding of FF parameters for bonds works
        self.top.add_ff_params()
        self.assertEqual(len(self.top.molecules[0].subsections[2].entries[3].params_state_a), 2)

    def test_add_params_angle(self):
        # checks if explicit adding of FF parameters for angles works
        self.top.add_ff_params()
        self.assertEqual(len(self.top.molecules[0].subsections[4].entries[3].params_state_a), 4)

    def test_energy(self):
        # checks if reading/writing does not change total energy
        ref_ener = gml.calc_gmx_energy('pentapeptide.pdb', 'pentapeptide.top')['potential']
        self.top.save_top('test.top')
        new_ener = gml.calc_gmx_energy('pentapeptide.pdb', 'test.top')['potential']
        os.remove('test.top')
        self.assertAlmostEqual(new_ener, ref_ener)

    def test_clear_sections(self):
        # checks if extra sections (like unused ions) are cleared correctly
        self.top.clear_sections()
        self.assertEqual(len(self.top.sections), 3)

    def test_system_properties(self):
        # checks if system properties are returned correctly (both values have to be 1)
        self.assertEqual(self.top.nmol('Protein_chain_A') * len(self.top.system), 1)

    def test_clear_params(self):
        # checks if parameters are cleared properly
        self.top.clear_ff_params()
        self.assertEqual(len(self.top.parameters.get_subsections('dihedraltypes')[0].entries_param), 95)

    def test_list_bonds(self):
        # checks if bonds are listed
        self.assertEqual(len(self.top.molecules[0].list_bonds(returning=True)), 87)

    def test_add_bond(self):
        # checks correct addition of bonds
        self.top.molecules[0].add_bond(1, 40)
        self.assertEqual(len(self.top.molecules[0].list_bonds(returning=True)), 88)

    def test_clone_type(self):
        # checks if cloning a type does not change energy
        ref_ener = gml.calc_gmx_energy('pentapeptide.pdb', 'pentapeptide.top')['potential']
        self.top.parameters.clone_type('H', new_type='HY')
        self.top.molecules[0].set_type('HY', atomname='HN')
        self.top.save_top('test.top')
        new_ener = gml.calc_gmx_energy('pentapeptide.pdb', 'test.top')['potential']
        os.remove('test.top')
        self.assertAlmostEqual(new_ener, ref_ener)


    def test_init_valid_file(self):
        self.assertIsInstance(self.top, gml.Top)

    def test_init_gmx_dir(self):
        self.assertEqual(self.top.gromacs_dir, '/usr/share/gromacs/top')

    def test_init_gmx_exe(self):
        self.assertEqual(self.top.gmx_exe, '/usr/bin/gmx')

    def test_init_pdb(self):
        self.top.add_pdb('pentapeptide.pdb')
        self.assertIsInstance(self.top.pdb, gml.Pdb)

#    def test_init_ignore_ifdef(self):
#        top = gml.Top('pentapeptide.top', ignore_ifdef=True)
#        self.assertTrue(self.top.ignore_ifdef)

    def test_init_define(self):
        top = gml.Top('pentapeptide.top', define={'TEST': 1})
        self.assertEqual(top.defines['TEST'], 1)

#    def test_init_keep_all(self):
#        top = gml.Top('pentapeptide.top', keep_all=True)
#        self.assertTrue(top.keep_all)

    def test_init_suppress(self):
        top = gml.Top('pentapeptide.top', suppress=True)
        self.assertTrue(top.suppress)

#    def test_init_amber(self):
#        top = gml.Top('', amber=True)
#        self.assertTrue('#define _FF_AMBER' in top._contents[0])
#
#    def test_init_charmm(self):
#        top = gml.Top('', charmm=True)
#        self.assertTrue('#define _FF_CHARMM' in top._contents[0])
#
#    def test_from_selection_valid(self):
#        selection = 'resid 1'
#        subtop = self.top.from_selection(selection)
#        self.assertIsInstance(subtop, gml.Top)
#
#    def test_from_selection_invalid(self):
#        selection = 'invalid selection'
#        with self.assertRaises(RuntimeError):
#            self.top.from_selection(selection)
#
#    def test_save_from_selection_valid(self):
#        selection = 'resid 1'
#        self.top.save_from_selection(selection, 'subsystem.top')
#        self.assertTrue(os.path.exists('subsystem.top'))
#        os.remove('subsystem.top')
#
#    def test_save_from_selection_invalid(self):
#        selection = 'invalid selection'
#        with self.assertRaises(RuntimeError):
#            self.top.save_from_selection(selection, 'subsystem.top')
#
#    def test_init_pdb(self):
#        top = gml.Top('pentapeptide.top', pdb='pentapeptide.pdb')
#        self.assertIsInstance(top.pdb, gml.Pdb)
#
#    def test_add_pdb_invalid(self):
#        with self.assertRaises(FileNotFoundError):
#            self.top.add_pdb('invalid.pdb')

#    def test_add_ff_params_valid(self):
#        self.top.add_ff_params('bonds')
#        self.assertIsEqual(len(self.top.molecules[0].bonds_section.entries_bonded[0].params_state_a), 2)

#    def test_add_ff_params_invalid(self):
#        with self.assertRaises(KeyError):
#            self.top.add_ff_params('invalid_section')



#    def test_load_frcmod_valid(self):
#        self.top.load_frcmod('valid.frcmod')
#        self.assertIsNotNone(self.top.parameters.get_subsection('bonds'))
#
#    def test_load_frcmod_invalid(self):
#        with self.assertRaises(FileNotFoundError):
#            self.top.load_frcmod('invalid.frcmod')
#
#    def test_add_molecule_from_file_valid(self):
#        self.top.add_molecule_from_file('valid.itp')
#        self.assertIsNotNone(self.top.get_molecule('NEW_MOL'))
#
#    def test_add_molecule_from_file_invalid(self):
#        with self.assertRaises(FileNotFoundError):
#            self.top.add_molecule_from_file('invalid.itp')
#
#    def test_add_molecule_from_file_molnames(self):
#        self.top.add_molecule_from_file('valid.itp', molnames=['NEW_MOL'])
#        self.assertIsNotNone(self.top.get_molecule('NEW_MOL'))
#
#    def test_add_molecule_from_file_molcount(self):
#        self.top.add_molecule_from_file('valid.itp', molcount=5)
#        self.assertEqual(self.top.nmol('NEW_MOL'), 5)
#
#    def test_add_molecule_from_file_prefix_type(self):
#        self.top.add_molecule_from_file('valid.itp', prefix_type='TEST')
#        mol = self.top.get_molecule('NEW_MOL')
#        self.assertTrue('TEST' in mol.atoms[0].type)
#
#    def test_add_parameters_from_file_valid(self):
#        self.top.add_parameters_from_file('valid.itp')
#        self.assertIsNotNone(self.top.parameters.get_subsection('bonds'))
#
#    def test_add_parameters_from_file_invalid(self):
#        with self.assertRaises(FileNotFoundError):
#            self.top.add_parameters_from_file('invalid.itp')
#
#    def test_add_parameters_from_file_sections(self):
#        self.top.add_parameters_from_file('valid.itp', sections=['bonds'])
#        self.assertIsNotNone(self.top.parameters.get_subsection('bonds'))
#
#    def test_add_parameters_from_file_overwrite(self):
#        self.top.add_parameters_from_file('valid.itp', overwrite=True)
#        self.assertIsNotNone(self.top.parameters.get_subsection('bonds'))
#
#    def test_add_parameters_from_file_prefix_type(self):
#        self.top.add_parameters_from_file('valid.itp', prefix_type='TEST')
#        bonds = self.top.parameters.get_subsection('bonds')
#        self.assertTrue('TEST' in bonds.entries_param[0].types[0])


if __name__ == "__main__":
    unittest.main()

