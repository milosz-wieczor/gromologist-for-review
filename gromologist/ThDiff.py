import os
import gromologist as gml
from typing import Union, Optional, TypeVar, Tuple
from itertools import combinations_with_replacement
from multiprocessing import Pool
from copy import deepcopy

gmlMod = TypeVar("gmlMod", bound="Mod")


class Mod:
    counter = 0  # only used for naming with unique numbers

    def __init__(self, top: Union[str, gml.Top], structure: str, selections: Optional[list] = None):
        """
        Base class to work with parameter modifications.
        The class is actually used for the calculation of original
        (unmodified) quantities, and is subclassed by ModParam,
        ModNbfix and ModAtom.
        :param top: list of Top objects corresponding to the modified topology
        :param structure: str, path to a structure file (.gro or .pdb) compatible with the topology
        :param selections: list of str, selections that will define atoms to be modified
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
        self.structure = structure
        self.topname = ''
        self.trajs = []
        self.dpar = 0.001
        self.sigma, self.chg, self.eps, self.nbs, self.nbe, self.dih, self.ang = [''] * 7
        self.atoms = [self.top.get_atoms(sel) for sel in selections]

    def __str__(self) -> str:
        """
        string representation of a Mod object, should be unique
        :return: mod one-letter code + Mod index, 1-based
        """
        return "{}-{}".format(''.join([self.sigma, self.chg, self.eps, self.nbs, self.nbe, self.dih, self.ang]),
                              self.counter)

    def __eq__(self, other: gmlMod) -> bool:
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
        final_name = self.topname + '/' + name + '-' + self.top.fname
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


class ModAtom(Mod):
    def __init__(self, top: Union[str, gml.Top], structure: str, selections: list, changes: str):
        """
        Subclass to work with modified atomic params (that is,
        sigma and epsilon defined for individual types, or charge defined
        for individual atoms).
        The class takes in a set of atoms that all correspond
        to a single atom type (has to be the same for all atoms);
        these atoms' type will be cloned from "type" to "Ytype",
        and the corresponding sigma/epsilon/charge parameters will be modified.
        :param top: a Top object corresponding to the modified topology
        :param selections: list of str, selections that will define atoms to be modified
        :param changes: single character denoting the change to be performed
        """
        super(ModAtom, self).__init__(top, structure, selections)
        self.sigma = 's' if 's' in changes else ''
        self.eps = 'e' if 'e' in changes else ''
        self.chg = 'c' if 'c' in changes else ''
        assert len(self.atoms) == 1
        self.atoms = self.atoms[0]
        if not (self.chg or self.eps or self.sigma):
            raise ValueError("no mode selected appropriate for an ModAtom object")
        self.types_are_consistent()
        self.type = self.atoms[0].type
        if not self.chg:
            self.top.parameters.clone_type(atomtype=self.type, prefix='Y')
        self.mod_sigma_eps()
        self.mod_chg()

    def mod_chg(self):
        """
        looks for a line in topology section 'atoms'
        and modifies it as needed
        :return: None
        """
        for atom in self.atoms:
            if "TD-MODDED" not in atom.comment:
                atom.comment = atom.comment + "TD-MODDED"
                atom.mass_b = atom.mass
                if self.chg:
                    atom.type_b = atom.type
                    atom.charge_b = atom.charge + self.dpar
                if self.sigma or self.eps:
                    atom.type_b = 'Y' + atom.type
                    atom.charge_b = atom.charge

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
        consistent = all([a.type == self.atoms[0].type for a in self.atoms])
        if not consistent:
            raise ValueError(
                "atoms within a single charge, sigma or epsilon modification need to have consistent types. "
                "To change atoms with different types, define multiple modifications")


class ModNbfix(Mod):
    def __init__(self, top: Union[str, gml.Top], structure: str, selections: list, changes: str):
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
        :param selections: list of str, selections that will define atoms to be modified
        :param changes: single character denoting the change to be performed
        """
        super(ModNbfix, self).__init__(top, structure, selections)
        self.nbs = 'n' if 'n' in changes else ''
        self.nbe = 'm' if 'm' in changes else ''
        if not (self.nbs or self.nbe):
            raise ValueError("no mode selected appropriate for an ModParam object")
        self.types = (self.atoms[0][0].type, self.atoms[1][0].type)
        if self.types[0] == self.types[1]:
            self.prefixes = ['Y', 'Y']
        else:
            self.prefixes = ['Q', 'Y']
        self.types_are_consistent()
        for type, prefix in zip(self.types, self.prefixes):
            if prefix + type not in self.top.defined_atomtypes:
                self.top.parameters.clone_type(atomtype=type, prefix=prefix)
        # we mod the first type
        self.prefix_type(self.atoms[0], self.types[0], self.prefixes[0])
        # and then second, if different than first
        if self.types[0] != self.types[1]:
            self.prefix_type(self.atoms[1], self.types[1], self.prefixes[1])
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

    def prefix_type(self, atoms: list, atom_type: str, prefix: str):
        """
        looks for a line in topology section 'atoms'
        and modifies it as needed
        :return: None
        """
        for atom in atoms:
            if "TD-MODDED" not in atom.comment:
                atom.mass_b = atom.mass
                atom.type_b = prefix + atom_type
                atom.charge_b = atom.charge
                atom.comment = atom.comment + "TD-MODDED"

    def types_are_consistent(self):
        """
        if many atoms are passed in a single input line, need to make sure
        they are all of the same type; raises ValueError if this is not the case
        :return: None
        """
        for atomset in self.atoms:
            consistent = all([atomset[0].type == atom.type for atom in atomset])
            if not consistent:
                raise ValueError("atoms within a single nbfix modification need to have consistent types")


class ModParam(Mod):
    def __init__(self, top: Union[str, gml.Top], structure: str, selections: list, changes: str,
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
        :param selections: list of str, selections that will define atoms to be modified
        :param changes: str, single character denoting the change to be performed
        """
        super(ModParam, self).__init__(top, structure, selections)
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
        self.types = [atoms[0].type for atoms in self.atoms]
        self.nums = [[atom.num for atom in atomset] for atomset in self.atoms]
        if self.ang:
            for molname in list({atom.molname for atomset in self.atoms for atom in atomset}):
                self.top.get_molecule(molname).add_ff_params('angles')
                self.mod_ang(molname)
        if self.dih:
            for molname in list({atom.molname for atomset in self.atoms for atom in atomset}):
                self.top.get_molecule(molname).add_ff_params('dihedrals')
                self.mod_dih(molname)

    def mod_ang(self, molname):
        """
        looks for lines in topology section 'angles'
        and modifies them as needed
        :return: None
        """
        subsection = self.top.get_molecule(molname).get_subsection('angles')
        assert len(self.types) == 3
        for entry in subsection.entries_bonded:
            entry.read_types()
            if all([i == j for i, j in zip(entry.types_state_a, self.types)]) or \
                    all([i == j for i, j in zip(entry.types_state_a, self.types[::-1])]):
                if all([i in numset for i, numset in zip(entry.atom_numbers, self.nums)] or \
                       [i in numset for i, numset in zip(entry.atom_numbers, self.nums[::-1])]):
                    if "TD-MODDED" not in entry.comment:
                        entry.params_state_b = entry.params_state_a[:]
                        entry.params_state_b[0] += self.dpar
                        entry.comment = entry.comment + "TD-MODDED"

    def mod_dih(self, molname):
        """
        looks for lines in topology section 'dihedrals'
        and modifies them as needed
        :return: None
        """
        subsection = self.top.get_molecule(molname).get_subsection('dihedrals')
        assert len(self.types) == 4
        for entry in subsection.entries_bonded:
            entry.read_types()
            if all([i == j for i, j in zip(entry.types_state_a, self.types)]) or \
                    all([i == j for i, j in zip(entry.types_state_a, self.types[::-1])]):
                if all([i in numset for i, numset in zip(entry.atom_numbers, self.nums)] or \
                       [i in numset for i, numset in zip(entry.atom_numbers, self.nums[::-1])]):
                    if entry.params_state_a[-1] == self.period:
                        if "TD-MODDED" not in entry.comment:
                            entry.params_state_b = entry.params_state_a[:]
                            entry.params_state_b[1] += self.dpar
                            entry.comment = entry.comment + "TD-MODDED"

    def types_are_consistent(self):
        """
        checks whether all input atoms share the same type signature,
        (t1 t2 t3 t4) for dihedrals or (t1 t2 t3) for angles;
        raises ValueError if not the case
        :return: None
        """
        for atomset in self.atoms:
            ok = all([atom.type == atomset[0].type for atom in atomset])
            if not ok:
                raise ValueError("atoms within a single parameter modification need to have consistent types")


class ThermoDiff:
    def __init__(self, temperature=300):
        self.mods = []
        self.temperature = temperature
        self.trajs = []
        self.path = None
        self.derivatives = {}
        self.discrete_free_energy_derivatives = {}
        self.profile_free_energy_derivatives = {}
        self.discrete_observable_derivatives = {}
        self.profile_observable_derivatives = {}

    def add_mod(self, top: Union[str, gml.Top], structure: str, modtype: str, selections: list):
        """
        Adds a new modified topology for derivative calculation
        :param top: str or gml.Top, the topology compatible with the system
        :param structure: str, a path to the Pdb or Gro file
        :param modtype: str, type of the modification, one of the following: c, s, e, n, m, a, d
        :param selections: list of str, selections that will define atoms to be modified
        :return: None
        """
        if modtype in 'cse':
            self.mods.append(ModAtom(top, structure, selections, modtype))
        elif modtype in 'ad':
            self.mods.append(ModParam(top, structure, selections, modtype))
        elif modtype in 'nm':
            self.mods.append(ModNbfix(top, structure, selections, modtype))

    def add_all_nbfix_mods(self, top: Union[str, gml.Top], structure: str):
        topology = gml.Top(top) if isinstance(top, str) else top
        topology.clear_sections()
        topology.clear_ff_params()
        for type_pair in sorted(list(combinations_with_replacement(topology.defined_atomtypes, 2))):
            print(f"Adding modification: {type_pair[0]}, {type_pair[1]}")
            self.add_mod(deepcopy(topology), structure, modtype='n',
                         selections=[f'type {type_pair[0]}', f'type {type_pair[1]}'])
            self.add_mod(deepcopy(topology), structure, modtype='m',
                         selections=[f'type {type_pair[0]}', f'type {type_pair[1]}'])

    def add_traj(self, top: Union[str, gml.Top], traj: str, datasets: Optional[dict] = None,
                 weights: Optional[list] = None):
        """
        Adds a trajectory on which the derivatives will be calculated,
        including pre-calculated data and frame-wise weights
        :param top: str or gml.Top, the topology compatible with the trajectory
        :param traj: str, path to the trajectory
        :param datasets: dict, named datasets in the format alias: list of [float/int/str]
        :param weights: list of float, per-frame weights (e.g. from Umbrella Sampling, Metadynamics or other reweighing
        schemes)
        :return: None
        """
        if isinstance(top, str):
            top = gml.Top(top)
        datasets = datasets if datasets is not None else {}
        weights = weights if weights is not None else []
        self.trajs.append({'path': traj, 'top': top, 'id': len(self.trajs), 'datasets': datasets, 'weights': weights})

    def get_traj(self, trajid: int) -> dict:
        """
        Gets a trajectory by its integer ID
        :param trajid: int, ID of the trajectory
        :return: dict, the trajectory and its associated data
        """
        return [traj for traj in self.trajs if traj['id'] == trajid][0]

    def add_data_to_traj(self, trajid: int, datasets: dict):
        """
        If not specified in self.add_traj, data can be added to the trajectory
        in a separate step
        :param trajid: int, ID of the trajectory
        :param datasets: dict, named datasets in the format alias: list of [float/int/str]
        :return: None
        """
        traj = self.get_traj(trajid)
        traj['datasets'].update(datasets)

    def add_weights_to_traj(self, trajid: int, weights: list):
        """
        If not specified in self.add_traj, weights can be added to the trajectory
        in a separate step
        :param trajid: int, ID of the trajectory
        :param weights: list of float, per-frame weights (e.g. from Umbrella Sampling, Metadynamics or other reweighing
        schemes)
        :return: None
        """
        traj = self.get_traj(trajid)
        traj['weights'] = weights

    @staticmethod
    def equal_tops(top1: gml.Top, top2: gml.Top) -> bool:
        """
        Checks (in a quick and approximate way) if two topologies are identical
        :param top1: gml.Top, first top to check
        :param top2: gml.Top, second top to check
        :return: bool, whether top1 and top2 are identical
        """
        if all(a1.atomname == a2.atomname for a1, a2 in zip(top1.atoms, top2.atoms)):
            return True
        else:
            return False

    def prep_files(self):
        """
        Prepares the directory tree for the calculations
        :return:
        """
        try:
            os.mkdir('working')
        except FileExistsError:
            print("The 'working' directory already exists, will overwrite its contents; clean incomplete 'rerun.xvg' "
                  "files to avoid problems")
        self.path = os.getcwd()
        for mod in self.mods:
            try:
                os.mkdir(f'working/{str(mod)}')
            except FileExistsError:
                pass
            mod.save_mod(f'{self.path}/working/{str(mod)}', str(mod))

    def launch_reruns(self):
        """
        Launches reruns for each legal mod-trajectory combination
        :return: None
        """
        pairs = []
        for mod in self.mods:
            for traj in self.trajs:
                if self.equal_tops(mod.top, traj['top']):
                    datasets = self.get_traj(traj['id'])['datasets']
                    print(f"Calculating rerun for {str(mod)}, traj {traj['path']}")
                    pairs.append((mod, traj, datasets))
        with Pool() as p:
            deriv_dicts = p.map(self.launch_rerun, pairs)
        for dd in deriv_dicts:
            self.derivatives.update(dd)

    @staticmethod
    def launch_rerun(mod_traj_ds):
        """
        Launches an individual alchemical rerun and reads the data, can be parallelized
        (skips the calculation if already performed)
        :param mod_traj_ds: tuple, contains the mod, traj and dataset objects
        :return: None
        """
        derivatives = {}
        mod, traj, datasets = mod_traj_ds
        mod.goto_mydir()
        if 'rerun.xvg' not in os.listdir('.'):
            derivatives[(mod.counter, traj['id'])] = gml.calc_gmx_dhdl('../../' + mod.structure, str(mod) + '-' +
                                                                       mod.top.top, '../../' + traj['path'])
        else:
            derivatives[(mod.counter, traj['id'])] = gml.read_xvg('rerun.xvg', [0])
        for key in datasets.keys():
            if not len(derivatives[(mod.counter, traj['id'])]) == len(datasets[key]):
                raise RuntimeError(f"The number of points in dataset {key} for trajectory num {traj['id']} "
                                   f"({len(datasets[key])}) does not match the derivatives "
                                   f"({len(derivatives[(mod.counter, traj['id'])])})")
        return derivatives

    def calc_weights(self):
        """
        Normalizes weights (or adds constant weights if not specified)
        :return: None
        """
        for traj in self.trajs:
            if not traj['weights']:
                traj['weights'] = [1] * self.get_length(traj['id'])
            wsum = sum(traj['weights'])
            traj['weights'] = [w/wsum for w in traj['weights']]

    def run(self):
        """
        Runs the set of calculations and prepares everything
        for data extraction
        :return: None
        """
        self.prep_files()
        self.launch_reruns()
        self.calc_weights()

    def get_mod(self, counter: int):
        """
        Gets the desired gml.Mod subclass instance by counter (ID)
        :param counter: int, ID of the modified topology
        :return: gml.Mod, selected modified topology
        """
        return [mod for mod in self.mods if mod.counter == counter][0]

    def get_length(self, trajid):
        """
        Gets the length of the trajectory from the length
        of the derivatives file (provided reruns have been performed)
        :param trajid: int, ID of the trajectory
        :return: int, length of the trajectory
        """
        for key in self.derivatives.keys():
            if key[1] == trajid:
                return len(self.derivatives[key])

    def dataset_is_numeric(self, dataset: str):
        """
        Checks if the user-specified dataset is numeric or string-based
        (in the latter case discrete states will be used)
        :param dataset: str, alias of the dataset
        :return: bool, whether the dataset is numeric
        """
        flat_dataset = {data for traj in self.trajs for data in traj['datasets'][dataset]}
        if all([isinstance(x, str) for x in flat_dataset]):
            return False
        elif all([isinstance(x, int) or isinstance(x, float) for x in flat_dataset]):
            return True
        else:
            raise RuntimeError(f"Dataset {dataset} contains mixed numeric and non-numeric entries")

    def get_flat_data(self, binning_dataset: str, deriv_dataset: str, mod) -> Tuple[list, list, list, list]:
        """
        Selects all numeric datasets that correspond to a given alias
        and provides flat (1-D) lists for data, weights, and derivatives
        :param binning_dataset: str, alias of the dataset used for binning/thresholding
        :param deriv_dataset: str, alias of the dataset used for derivative calculation (can be same as binning_dataset)
        :param mod: gml.Mod
        :return: tuple of lists: data, weights, and derivatives
        """
        counters_trajids = list(self.derivatives.keys())
        active_dataset_trajids = [traj['id'] for traj in self.trajs if binning_dataset in traj['datasets'].keys()
                                  and deriv_dataset in traj['datasets'].keys()]
        active_counters_trajids = [ct for ct in counters_trajids if ct[1] in active_dataset_trajids
                                   and ct[0] == mod.counter]
        flat_bin_data = [data for ct in active_counters_trajids for data in
                         self.get_traj(ct[1])['datasets'][binning_dataset]]
        flat_deriv_data = [data for ct in active_counters_trajids for data in
                           self.get_traj(ct[1])['datasets'][deriv_dataset]]
        flat_weights = [wght for ct in active_counters_trajids for wght in self.get_traj(ct[1])['weights']]
        flat_derivs = [der for ct in active_counters_trajids for der in self.derivatives[ct]]
        return flat_bin_data, flat_deriv_data, flat_weights, flat_derivs

    def calc_discrete_derivatives(self, dataset: Optional[str] = None, free_energy: Optional[bool] = True,
                                  threshold: Optional[list] = None, cv_dataset: Optional[str] = None):
        """
        Calculates selected derivatives for a number of discrete states
        :param dataset: str, alias for the dataset
        :param free_energy: bool, whether the derivative should be of free energy or of observable/CV
        :param threshold: list, for continues dataset specifies 2N boundaries used to define N discrete states
        :param cv_dataset: str, if specified then derivatives will be calculated based on 'dataset' but binning will be
        performed based on 'cv_dataset'
        :return: None
        """
        if dataset is None:
            print("No dataset label chosen, will assume each trajectory defines a separate state")
            for traj in self.trajs:
                traj['datasets']['discrete'] = [str(traj['id']) for _ in range(self.get_length(traj['id']))]
            dataset = 'discrete'
        if self.dataset_is_numeric(dataset) and threshold is None:
            raise RuntimeError("When calculating a derivative with numeric datasets, specify a threshold; when "
                               "calculating a derivative over discrete states, use strings as state labels")
        if threshold is not None:
            if len(threshold) % 2 == 1:
                raise RuntimeError("The list 'threshold' has to have an even number of entries, one for each starting "
                                   "and ending point of each state")
        states = [(round(x, 6), round(y, 6)) for x, y in zip(threshold[::2], threshold[1::2])] if threshold is not None else \
            sorted(list({data for traj in self.trajs for data in traj['datasets'][dataset]}))
        binning_dset, deriv_dset = (dataset, dataset) if cv_dataset is None else (cv_dataset, dataset)
        for mod in self.mods:
            binning_data, deriv_data, weights, derivs = self.get_flat_data(binning_dset, deriv_dset, mod)
            mean_derivatives = {st: 0 for st in states + [None]}
            mean_product = {st: 0 for st in states + [None]}
            mean_data = {st: 0 for st in states + [None]}
            for n in range(len(binning_data)):
                state_index = None
                if threshold is not None:
                    for x, y in zip(threshold[::2], threshold[1::2]):
                        if x <= binning_data[n] < y:
                            state_index = (round(x, 6), round(y, 6))
                            break
                else:
                    state_index = binning_data[n]
                mean_derivatives[state_index] += weights[n] * derivs[n]
                if not free_energy:
                    if not free_energy:
                        mean_product[state_index] += deriv_data[n] * weights[n] * derivs[n]
                        mean_data[state_index] += deriv_data[n]
            if free_energy:
                self.discrete_free_energy_derivatives[(str(mod), dataset)] = [mean_derivatives[x] / mod.dpar
                                                                              for x in mean_derivatives.keys()]
            else:
                mean_obs = {key: 0 for key in mean_derivatives.keys()}
                for key in mean_derivatives.keys():
                    mean_obs[key] = (1/0.008314*self.temperature) * (
                            mean_data[key] * mean_derivatives[key] - mean_product[key])
                self.discrete_observable_derivatives[(str(mod), dataset)] = [mean_obs[x] / mod.dpar
                                                                             for x in mean_obs.keys()]

    def calc_profile_derivatives(self, dataset: str, free_energy: Optional[bool] = True, nbins: Optional[int] = 50,
                                 cv_dataset: Optional[str] = None):
        """
        Calculates selected derivatives of a profile
        :param dataset: str, alias for the dataset
        :param free_energy: bool, whether the derivative should be of free energy or of observable/CV
        :param nbins: int, number of bins that will be used to calculate the profile
        :param cv_dataset: str, if specified then derivatives will be calculated based on 'dataset' but binning will be
        performed based on 'cv_dataset'
        :return: None
        """
        binning_dset, deriv_dset = (dataset, dataset) if cv_dataset is None else (cv_dataset, dataset)
        for mod in self.mods:
            binning_data, deriv_data, weights, derivs = self.get_flat_data(binning_dset, deriv_dset, mod)
            dmin, dmax = min(binning_data), max(binning_data)
            thresh = [dmin + n * (dmax-dmin)/(nbins+1) for n in range(nbins+1)]
            mean_derivatives = {0.5 * (x + y): 0 for x, y in zip(thresh[:-1], thresh[1:])}
            mean_product = {0.5 * (x + y): 0 for x, y in zip(thresh[:-1], thresh[1:])}
            mean_data = {0.5 * (x + y): 0 for x, y in zip(thresh[:-1], thresh[1:])}
            for n in range(len(binning_data)):
                for x, y in zip(thresh[:-1], thresh[1:]):
                    if x <= binning_data[n] < y:
                        mean_derivatives[0.5*(x+y)] += weights[n] * derivs[n]
                        if not free_energy:
                            mean_product[0.5*(x+y)] += deriv_data[n] * weights[n] * derivs[n]
                            mean_data[0.5*(x+y)] += deriv_data[n]
            if free_energy:
                self.profile_free_energy_derivatives[(str(mod), dataset)] = ([(x+y)/2 for x, y in zip(thresh[:-1],
                                                                                                      thresh[1:])],
                                                                             [mean_derivatives[x] / mod.dpar
                                                                             for x in mean_derivatives.keys()])
            else:
                mean_obs = {key: 0 for key in mean_derivatives.keys()}
                for key in mean_derivatives.keys():
                    mean_obs[key] = (1 / 0.008314 * self.temperature) * (
                                mean_data[key] * mean_derivatives[key] - mean_product[key])
                self.profile_observable_derivatives[(str(mod), dataset)] = ([(x+y)/2 for x, y in zip(thresh[:-1],
                                                                                                     thresh[1:])],
                                                                            [mean_obs[x] / mod.dpar
                                                                            for x in mean_obs.keys()])
