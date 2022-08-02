import os
import gromologist as gml
from typing import Union, Optional, Type, Iterable


class Mod:
    counter = 0  # only used for naming with unique numbers

    def __init__(self, top: Union[str, gml.Top], resnr: Optional[list] = None, resnames: Optional[list] = None,
                 names: Optional[list] = None, molname: Optional[str] = None):
        """
        Base class to work with parameter modifications.
        The class is actually used for the calculation of original
        (unmodified) quantities, and is subclassed by ModParam,
        ModNbfix and ModAtom.
        :param top: list of Top objects corresponding to the modified topology
        :param resnr: a list of tuples, each containing a list of 3 or 4 residue numbers
        :param resnames: a list of tuples, each containing a list of 3 or 4 residue names
        :param names: a list of tuples, each containing a list of 3 or 4 atom names
        """
        Mod.counter += 1
        self.counter = Mod.counter
        if isinstance(top, str):
            self.top = gml.Top(top)
        elif isinstance(top, gml.Top):
            self.top = top
        else:
            raise RuntimeError('"top" has to be either a path to file or a gml.Top object')
        self.data = None
        self.diff = None
        self.molname = molname
        self.topname = ''
        self.trajs = []
        self.dpar = 0.001
        self.sigma, self.chg, self.eps, self.nbs, self.nbe, self.dih, self.ang, self.nothing = [''] * 8
        assert bool(names) == bool(resnames) == bool(resnr)
        self.names = names
        self.res = resnr
        self.resnames = resnames
        assert len(self.names) == len(self.res) == len(self.resnames) > 0

    def __str__(self) -> str:
        """
        string representation of a Mod object, should be unique
        :return: mod one-letter code + Mod index, 1-based
        """
        return "{}-{}".format(''.join([self.sigma, self.chg, self.eps, self.nbs, self.nbe, self.dih, self.ang,
                                       self.nothing]), self.counter)

    def __eq__(self, other: Type[gml.Mod]) -> bool:
        """
        assumes unique numbering of instances
        :param other: another Mod (or subclass) instance to compare with
        :return: bool
        """
        return str(self) == str(other)

    def save_mod(self, path: str, name: str):
        """
        wraps Top's save_mod to catch the file name
        :param path: str, filename to save .top to
        :param name: str, name of the file
        :return: None
        """
        self.topname = path
        final_name = self.topname + name + '-' + self.top.fname
        if not os.path.exists(final_name):
            self.top.save_top(final_name)
        else:
            print(f'file {final_name} already exists, leaving as it is')

    def goto_mydir(self):
        """
        shortcut to the ModAtom's working directory
        :return: None
        """
        os.chdir(self.topname)

    def get_type(self, name: str, res: int, resname: str, molname: str) -> str:
        mol = self.top.get_molecule(molname)
        return mol.get_atom(f'name {name} and resid {res} and resname {resname}').type


class ModAtom(Mod):
    def __init__(self, top: Union[str, gml.Top], resnr: list, resnames: list, names: list, molname: str, changes: str):
        """
        Subclass to work with modified atomic params (that is,
        sigma and epsilon defined for individual types, or charge defined
        for individual atoms).
        The class takes in a set of atoms that all correspond
        to a single atom type (has to be the same for all atoms);
        these atoms' type will be cloned from "type" to "Ytype",
        and the corresponding sigma/epsilon/charge parameters will be modified.
        :param top: a Top object corresponding to the modified topology
        :param resnr: a list of residue numbers
        :param resnames: a list of residue names
        :param names: a list of of atom names
        :param changes: single character denoting the change to be performed
        """
        super(ModAtom, self).__init__(top, resnr, resnames, names, molname)
        self.sigma = 's' if 's' in changes else ''
        self.eps = 'e' if 'e' in changes else ''
        self.chg = 'c' if 'c' in changes else ''
        if not (self.chg or self.eps or self.sigma):
            raise ValueError("no mode selected appropriate for an ModAtom object")
        if names and resnr:
            self.types_are_consistent()
            self.type = self.get_type(self.names[0], self.res[0], self.resnames[0], self.molname)
            if not self.chg:
                self.top.parameters.clone_type(atomtype=self.type, prefix='Y')
            self.mod_sigma_eps()
            self.mod_chg()
        else:
            self.names = ['orig']
            self.res = []

    def mod_chg(self):
        """
        looks for a line in topology section 'atoms'
        and modifies it as needed
        :return: None
        """
        for x in range(len(self.names)):
            name = self.names[x]
            resname = self.resnames[x]
            resnum = self.res[x]
            subsection = self.top.get_molecule(self.molname).get_subsection('atoms')
            for e in subsection.entries_atom:
                if e.atomname == name and e.resid == resnum and e.resname == resname:
                    e.mass_b = e.mass
                    if self.chg:
                        e.type_b = e.type
                        e.charge_b = e.charge + self.dpar
                    if self.sigma or self.eps:
                        e.type_b = 'Y' + e.type
                        e.charge_b = e.charge

    def mod_sigma_eps(self):
        """
        looks for a line in topology section 'atomtypes'
        and modifies it as needed
        :return: None
        """
        subsection = self.top.parameters.get_subsection('atomtypes')
        for e in subsection.entries_param:
            if e.types[0] == 'Y' + self.type:
                if self.sigma:
                    e.params[0] += self.dpar
                if self.eps:
                    e.params[1] += self.dpar

    def types_are_consistent(self):
        """
        if many atoms are passed in a single input line, need to make sure
        they are all of the same type
        :return: True if all are of the same type, False otherwise
        """
        consistent = all([self.get_type(self.names[0], self.res[0], self.resnames[0], self.molname)
                          == self.get_type(self.names[n], self.res[n], self.resnames[n], self.molname)
                          for n in range(len(self.names))])
        if not consistent:
            raise ValueError(
                "atoms within a single charge, sigma or epsilon modification need to have consistent types")


class ModNbfix(Mod):
    def __init__(self, top: Union[str, gml.Top], resnr: list, resnames: list, names: list, molname: str, changes: str):
        """
        Subclass to work with modified NBFIX params (pairwise corrections,
        sigma and epsilon defined for pairs of types).
        The class allows two lists (actually, tuples) of atoms that correspond
        to two atom types (atomtypes have to be consistent within each tuple);
        these atoms' types will be cloned from (type1, type2) to (Xtype1,
        Ytype2), and the corresponding NBFIX parameters will be modified.
        The modification can be applied to all atoms bearing the specific type
        as well as only a subgroup of these.
        :param top: a Top object corresponding to the modified topology
        :param resnr: a list of two tuples, each containing a list of residue numbers
        :param resnames: a list of two tuples, each containing a list of residue names
        :param names: a list of two tuples, each containing a list of atom names
        :param changes: single character denoting the change to be performed
        """
        super(ModNbfix, self).__init__(top, resnr, resnames, names, molname)
        self.nbs = 'n' if 'n' in changes else ''
        self.nbe = 'm' if 'm' in changes else ''
        if not (self.nbs or self.nbe):
            raise ValueError("no mode selected appropriate for an ModParam object")
        for tnames, tres, tresnames in zip(self.names, self.res, self.resnames):
            self.types_are_consistent(tnames, tres, tresnames)
            try:
                self.types = (self.get_type(self.top, tnames[0][0], tres[0][0], tresnames[0][0]),
                              self.get_type(self.top, tnames[1][0], tres[1][0], tresnames[1][0]))
            except IndexError:
                pass
            else:
                if self.types[0] == self.types[1]:
                    self.prefixes = ['Y', 'Y']
                else:
                    self.prefixes = ['Q', 'Y']
                # we mod the first type
                self.top.parameters.clone_type(atomtype=self.types[0], prefix=self.prefixes[0])
                self.prefix_type(tnames, tres, tresnames, 0)
                # and then second, if different than first
                if self.types[0] != self.types[1]:
                    self.top.parameters.clone_type(atomtype=self.types[1], prefix=self.prefixes[1])
                    self.prefix_type(tnames, tres, tresnames, 1)
                self.mod_nb_sigma_eps()

    def mod_nb_sigma_eps(self):
        """
        looks for a line in topology section 'nonbond_params'
        and modifies it as needed
        :return: None
        """
        subsection = self.top.parameters
        if self.nbs:
            subsection.add_nbfix(type1=self.prefixes[0] + self.types[0], type2=self.prefixes[1] + self.types[1],
                                 mod_sigma=self.dpar, action_default='m')
        elif self.nbe:
            subsection.add_nbfix(type1=self.prefixes[0] + self.types[0], type2=self.prefixes[1] + self.types[1],
                                 mod_epsilon=self.dpar, action_default='m')

    def prefix_type(self, tnames: Iterable, tres: Iterable, tresnames: Iterable, prefix_number: int):
        """
        looks for a line in topology section 'atoms'
        and modifies it as needed
        :return: None
        """
        for x in range(len(tnames[prefix_number])):
            name = tnames[prefix_number][x]
            resname = tresnames[prefix_number][x]
            resnum = tres[prefix_number][x]
            subsection = self.top.get_molecule(self.molname)
            subsection.set_type(resname=resname, resid=resnum, atomname=name, prefix=self.prefixes[prefix_number])

    def types_are_consistent(self, tnames: Iterable, tres: Iterable, tresnames: Iterable):
        """
        if many atoms are passed in a single input line, need to make sure
        they are all of the same type; raises ValueError if this is not the case
        :return: None
        """
        consistent1 = all([self.get_type(tnames[0][0], tres[0][0], tresnames[0][0], self.molname)
                           == self.get_type(tnames[0][n], tres[0][n], tresnames[0][n], self.molname)
                           for n in range(len(tnames[0]))])
        consistent2 = all([self.get_type(tnames[1][0], tres[1][0], tresnames[1][0], self.molname)
                           == self.get_type(tnames[1][n], tres[1][n], tresnames[1][n], self.molname)
                           for n in range(len(tnames[1]))])
        if not consistent1 or not consistent2:
            raise ValueError("atoms within a single nbfix modification need to have consistent types")


class ModParam(Mod):
    def __init__(self, top: Union[str, gml.Top], resnr: list, resnames: list, names: list, molname: str, changes: str,
                 period: Optional[int] = -1):
        """
        Subclass to work with modified bonded params (angles and dihedrals).
        The class allows a lists of tuples of atoms that each corresponds
        to a single bonded term, angle or dihedral (all need to have
        the same sets of atomtypes for the calculation to be meaningful);
        the corresponding parameters will be modified in the molecule definition.
        The modification can be applied to all atoms bearing the specific type
        as well as only a subgroup of these.
        :param top: Top, a Top object corresponding to the modified topology
        :param resnr: list, a list of tuples, each containing a list of 3 or 4 residue numbers
        :param resnames: list, a list of tuples, each containing a list of 3 or 4 residue names
        :param names: list, a list of tuples, each containing a list of 3 or 4 atom names
        :param changes: str, single character denoting the change to be performed
        """
        super(ModParam, self).__init__(top, resnr, resnames, names, molname)
        self.dih = 'd' if 'd' in changes else ''  # implement add to top line
        self.ang = 'a' if 'a' in changes else ''  # implement add to top line
        self.period = period
        if self.dih and self.ang:
            raise RuntimeError("Please use separate Mods entries for angles and dihedrals")
        if not (self.dih or self.ang):
            raise ValueError("no mode selected appropriate for an ModParam object")
        if period < 0 and self.dih:
            raise RuntimeError("Please choose periodicity for the dihedral to be modified")
        self.types_are_consistent()
        self.types = tuple(self.get_type(self.top, n, r, rn) for n, r, rn in
                           zip(self.names[0], self.res[0], self.resnames[0]))
        if self.ang:
            self.top.get_molecule(self.molname).add_ff_params('angles')
            self.mod_ang()
        if self.dih:
            self.top.get_molecule(self.molname).add_ff_params('dihedrals')
            self.mod_dih()

    def mod_ang(self):
        """
        looks for lines in topology section 'angles'
        and modifies them as needed
        :return: None
        """
        subsection = self.top.get_molecule(self.molname).get_subsection('angles')
        assert len(self.types) == 3
        for entry in subsection.entries_bonded:
            entry.read_types()
            if all([i == j for i, j in zip(entry.types_state_a, self.types)]) or\
                    all([i == j for i, j in zip(entry.types_state_a, self.types[::-1])]):
                entry.params_state_b = entry.params_state_a[:]
                entry.params_state_b[0] += self.dpar

    def mod_dih(self):
        """
        looks for lines in topology section 'dihedrals'
        and modifies them as needed
        :return: None
        """
        subsection = self.top.get_molecule(self.molname).get_subsection('dihedrals')
        assert len(self.types) == 4
        for entry in subsection.entries_bonded:
            entry.read_types()
            if all([i == j for i, j in zip(entry.types_state_a, self.types)]) or \
                    all([i == j for i, j in zip(entry.types_state_a, self.types[::-1])]):
                if entry.params_state_a[-1] == self.period:
                    entry.params_state_b = entry.params_state_a[:]
                    entry.params_state_b[1] += self.dpar

    def types_are_consistent(self):
        """
        checks whether all input atoms share the same type signature,
        (t1 t2 t3 t4) for dihedrals or (t1 t2 t3) for angles;
        raises ValueError if not the case
        :return: None
        """
        types_list = [tuple(self.get_type(qn, qrn, qr, self.molname) for qn, qrn, qr in zip(n, rn, r))
                      for n, rn, r in zip(self.names, self.res, self.resnames)]
        ok = all([types_list[0] == types_list[n] for n in range(1, len(types_list))])
        if not ok:
            raise ValueError("atoms within a single parameter modification need to have consistent types")


class ThermoDiff:
    def __init__(self):
        self.mods = []
        self.trajs = []

    def add_mod(self, top: Union[str, gml.Top], modtype: str, resids: list, resnames: list, atomnames: list, molname: str):
        if modtype in 'cse':
            self.mods.append(ModAtom(top, resids, resnames, atomnames, molname, modtype))
        elif modtype in 'ad':
            self.mods.append(ModParam(top, resids, resnames, atomnames, molname, modtype))
        elif modtype in 'nm':
            self.mods.append(ModNbfix(top, resids, resnames, atomnames, molname, modtype))

    def add_traj(self, top: Union[str, gml.Top], traj: str):
        if isinstance(top, str):
            top = gml.Top(top)
        self.trajs.append((top, traj))

    @staticmethod
    def equal_tops(top1: gml.Top, top2: gml.Top) -> bool:
        if all(a1.atomname == a2.atomname for a1, a2 in zip(top1.atoms, top2.atoms)):
            return True
        else:
            return False
