import gromologist as gml


class Subsection:
    counter = {}
    
    def __init__(self, content, section):
        """
        Here we want to have:
          - a unique representation of the section (header/ID)
          - a list of entries
          - a binding to the Top class that holds the Section object

        :param content: list of strings, entire content of the section
        :param section: a Section instance that contains this Subsection
        """
        self.section = section
        self.conditional = False
        self.header = content[0].strip().strip('[]').strip()
        if ';' in self.header:
            pos = self.header.index(';')
            self.header = self.header[:pos].strip().strip('[]').strip()
        if self.header in Subsection.counter.keys():
            Subsection.counter[self.header] += 1
        else:
            Subsection.counter[self.header] = 1
        self.id = Subsection.counter[self.header]
        self.entries = []
        for element in content:
            if issubclass(type(element), gml.Entry):
                self.entries.append(element)
            elif isinstance(element, str) and element.strip() and not element.strip().startswith('['):
                self.entries.append(self.yield_entry(element))
        
    def yield_entry(self, line):
        """
        Decides which Entry subclass to return
        based on which Subsection subclass evokes this fn
        :param line: str, a line to be converted into an Entry instance
        :return: Entry, an instance of the proper Entry subclass
        """
        if line.strip()[0] in [';', '#']:
            return gml.Entry(line, self)
        elif isinstance(self, SubsectionParam):
            return gml.EntryParam(line, self)
        elif isinstance(self, SubsectionBonded):
            return gml.EntryBonded(line, self)
        elif isinstance(self, SubsectionAtom):
            return gml.EntryAtom(line, self)
        elif isinstance(self, Subsection):
            return gml.Entry(line, self)
    
    def __str__(self):
        """
        As section headers can be repeated, each section is denoted
        by a header and ID (ID corresponding to the consecutive numbering
        of the specific header, e.g. bonds-3 is the third "bonds" section
        :return: str, section label
        """
        return "{}-{}".format(self.header, self.id)
    
    def __repr__(self):
        return "Subsection {}".format(self.header, self.id)
    
    def __len__(self):
        return len(self.entries)
    
    def __iter__(self):
        """
        Useful if we want to iterate over entries as "for entry in subsection",
        allows us to mark self._entries as private
        :return: self
        """
        self.n = 0
        return self
    
    def __next__(self):
        n = self.n
        self.n += 1
        try:
            return self.entries[n]
        except IndexError:
            raise StopIteration

    def add_entry(self, new_entry, position=None):
        """
        Adds a single entry to the subsection, either at the end
        or in a specified position
        :param new_entry: gml.Entry, entry to be added
        :param position: where to add the entry (None is at the end)
        :return: None
        """
        new_entry.subsection = self
        if position is not None:
            position = int(position)
            self.entries.insert(position, new_entry)
        else:
            self.entries.append(new_entry)
    
    def add_entries(self, new_entries_list, position=None):
        """
        Adds multiple entries to the subsection, either at the end
        or in a specified position
        :param new_entries_list: list of gml.Entry, entries to be added
        :param position: where to add the entries (None is at the end)
        :return: None
        """
        for new_entry in new_entries_list:
            new_entry.subsection = self
        if position is not None:
            position = int(position)
            for new_entry in new_entries_list:
                self.entries.insert(position, new_entry)
                position += 1
        else:
            self.entries.extend(new_entries_list)
    
    def set_entry(self, line_number, new_line):
        """
        Sets content of a specified entry
        :param line_number: int, which entry to modify
        :param new_line: str, new content of the entry
        :return: None
        """
        self.entries[line_number] = new_line
    
    def get_entry(self, line_number):
        """
        Returns entry specified by line number
        :param line_number: int, which entry to return
        :return: str, subsection entry
        """
        return self.entries[line_number]

    def remove_entry(self, entry):
        self.entries.remove(entry)
        
        
class SubsectionBonded(Subsection):
    """
    SubsectionBonded contains a subsection with entries corresponding to bonded terms,
    e.g., bonds or dihedrals; should be included in SectionMol
    """
    n_atoms = {'bonds': 2, 'pairs': 2, 'angles': 3, 'dihedrals': 4,
               'cmap': 5, 'settles': 1, 'exclusions': 2, 'position_restraints': 1,
               'virtual_sites2': 3, 'constraints': 2}
    
    def __init__(self, content, section):
        super().__init__(content, section)
        self.bkp_entries = None
        self.atoms_per_entry = SubsectionBonded.n_atoms[self.header]
        self.prmtypes = self._check_parm_type()
        self.fstring = "{:5} " * (SubsectionBonded.n_atoms[self.header] + 1) + '\n'
    
    def __repr__(self):
        return "Subsection {} with interaction type {}".format(self.header, ' '.join(self.prmtypes))

    @property
    def label(self):
        return '{}-{}'.format(self.header, '_'.join(self.prmtypes))

    @property
    def entries_bonded(self):
        """
        Lists all entries that are of EntryBonded type (i.e. excluding comments,
        headers, empty lines etc)
        :return: list of EntryBonded instances
        """
        return [e for e in self.entries if isinstance(e, gml.EntryBonded)]
    
    def sort(self):
        """
        In case we want to sort entries after some are added at the end of the section
        :return: None
        """
        self.entries.sort(key=self._sorting_fn)
    
    def _sorting_fn(self, entry):
        """
        Comments should go first, then we sort based on first, second,
        ... column of the section
        :param entry: Entry, entry to be sorted
        :return: int, ordering number
        """
        if isinstance(entry, gml.Entry):
            return -1
        val = sum([i * 10**(4*(self.atoms_per_entry - n)) for n, i in enumerate(entry.atom_numbers)])
        return val

    def explicit_defines(self):
        """
        Substitutes predefined fields with their corresponding values
        as specified by the #DEFINE preprocessor commands
        :return: None
        """
        for entry in self.entries_bonded:
            entry.explicit_defines()
    
    def add_ff_params(self, force_all=False):
        """
        Adds explicit values of FF parameters in sections 'bonds', 'angles', 'dihedrals'
        :param force_all: bool, whether to add params from scratch if some are assigned already
        :return: None
        """
        matchings = {'bonds': 'bondtypes', 'angles': 'angletypes', 'dihedrals': 'dihedraltypes'}
        subsect_params = [sub for sub in self.section.top.parameters.subsections
                          if sub.header == matchings[self.header]]
        self.bkp_entries = self.entries[:]  # we can't change what we're iterating over, so we modify the copy
        for entry in self.entries_bonded:
            # let's omit the ones that already have params assigned:
            if entry.params_state_a:
                if not force_all:
                    continue
                else:
                    entry.params_state_a = []
                    entry.params_state_b = []
                    entry.params_state_a_entry = []
                    entry.params_state_b_entry = []
            self._add_ff_params_to_entry(entry, subsect_params)
        self.entries = self.bkp_entries[:]  # now restore the modified copy

    def find_used_ff_params(self):
        """
        Looks up a list of FF parameters that are actually used in the molecular
        system at hand
        :return: list of str, labels of used parameters
        """
        used_parm_entries = []
        matchings = {'bonds': 'bondtypes', 'angles': 'angletypes', 'dihedrals': 'dihedraltypes', 'cmap': 'cmaptypes'}
        subsect_params = [sub for sub in self.section.top.parameters.subsections if
                          sub.header == matchings[self.header]]
        for entry in self.entries_bonded:
            if not entry.params_state_a:
                used_parm_entries.extend(self._find_used_ff_params(entry, subsect_params))
        return used_parm_entries

    @staticmethod
    def _find_used_ff_params(entry, subsect_params):
        """
        Identifies the FF parameter(s) that is used for the particular
        bonded interaction entry
        :param entry: a BondedEntry instance, an entry we're finding the parameter for
        :param subsect_params: a SubsectionParam instance containing candidate parameters for the entry
        :return: list of str, labels of used parameters
        """
        entries = []
        int_type = entry.interaction_type
        entry.read_types()
        for types in [entry.types_state_a, entry.types_state_b]:
            wildcard_present = []
            non_wildcard_present = []
            for subsections in subsect_params:
                for parm_entry in [e for e in subsections.entries_param]:
                    if parm_entry.match(types, int_type):
                        is_wildcard = 'X' in parm_entry.types
                        if not wildcard_present and not is_wildcard:
                            entries.append(parm_entry.identifier)
                            non_wildcard_present += parm_entry.types
                        elif not wildcard_present and is_wildcard and not non_wildcard_present:
                            entries.append(parm_entry.identifier)
                            wildcard_present = parm_entry.types
                        elif wildcard_present and not is_wildcard:
                            raise RuntimeError("Wildcard ('X') entries were found prior to regular ones, please fix"
                                               "your FF parameters")
                        elif wildcard_present and is_wildcard:  # only add if multiple entries per given wildcard
                            if parm_entry.types == wildcard_present:
                                entries.append(parm_entry.identifier)
                            else:
                                pass
        return entries

    def find_missing_ff_params(self, fix_by_analogy=False, fix_B_from_A=False, fix_A_from_B=False, fix_dummy=False,
                               once=False):
        """
        Identifies FF parameters that cannot by matched to the existing set of
        bonded interactions, allowing to find & target these for fixing
        :param fix_by_analogy: dict, tells how to choose new parameters from existing equivalent types
        :param fix_B_from_A: bool, chooses missing parameters assuming B-state should have similar types to A-state
        :param fix_A_from_B: bool, chooses missing parameters assuming A-state should have similar types to B-state
        :param fix_dummy: bool, adds dummy (all-0) parameters
        :param once: bool, only shows the missing parameter once if repeated
        :return: None
        """
        matchings = {'bonds': 'bondtypes', 'angles': 'angletypes', 'dihedrals': 'dihedraltypes'}
        # TODO add cmaps?
        subsect_params = [sub for sub in self.section.top.parameters.subsections if
                          sub.header == matchings[self.header]]
        for entry in self.entries_bonded:
            self._find_missing_ff_params(entry, subsect_params, fix_by_analogy, fix_B_from_A, fix_A_from_B, fix_dummy,
                                         once)

    def _find_missing_ff_params(self, entry, subsect_params, fix_by_analogy, fix_B_from_A, fix_A_from_B, fix_dummy,
                                once):
        if (entry.params_state_a and entry.params_state_b) or (entry.params_state_a and not entry.types_state_b):
            return
        int_type = entry.interaction_type
        entry.read_types()
        found_a, found_b = False, False
        for subsections in subsect_params:
            for parm_entry in subsections.entries_param:
                if parm_entry.match(entry.types_state_a, int_type):
                    found_a = True
        if not found_a and not entry.params_state_a:
            if not once or entry.types_state_a not in self.section.printed:
                print(f'Couldn\'t find params for interaction type {entry.subsection.header} {int_type}, '
                      f'atom types {entry.types_state_a}, atom numbers {entry.atom_numbers}')
            if once and entry.types_state_a not in self.section.printed:
                self.section.printed.append(entry.types_state_a)
            if fix_by_analogy:
                candid = self._fix_by_analogy(fix_by_analogy, entry.types_state_a, subsect_params, int_type)
                if candid:
                    entry.params_state_a = candid
            if fix_dummy:
                entry.params_state_a = [0.0 for _ in gml.EntryBonded.fstr_suff[(self.header, self.prmtypes[0])]]
                print("setting dummy parameters: " + str(entry))
        if entry.types_state_b and not entry.params_state_b:
            for subsections in subsect_params:
                for parm_entry in subsections.entries_param:
                    if parm_entry.match(entry.types_state_b, int_type):
                        found_b = True
            if not found_b and not entry.params_state_b:
                if not once or entry.types_state_b not in self.section.printed:
                    print(f'Couldn\'t find params for interaction type {entry.subsection.header} {int_type}, '
                          f'atom types {entry.types_state_b}, atom numbers {entry.atom_numbers}')
                if once and entry.types_state_b not in self.section.printed:
                    self.section.printed.append(entry.types_state_b)
                if fix_by_analogy:
                    candid = self._fix_by_analogy(fix_by_analogy, entry.types_state_b, subsect_params, int_type)
                    if candid:
                        entry.params_state_b = candid
        if entry.types_state_b:
            if fix_B_from_A and not found_b and not entry.params_state_b:
                candid = self._fix_by_analogy({}, entry.types_state_b, subsect_params, int_type,
                                              other_typelist=entry.types_state_a)
                if candid:
                    entry.params_state_b = candid
            if fix_A_from_B and not found_a and not entry.params_state_a:
                candid = self._fix_by_analogy({}, entry.types_state_a, subsect_params, int_type,
                                              other_typelist=entry.types_state_b)
                if candid:
                    entry.params_state_a = candid

    @staticmethod
    def _fix_by_analogy(subst, typelist, subsect_params, int_type, other_typelist=None):
        new_params = []
        from_wildtype = None
        if other_typelist is None:
            types = [tp if tp not in subst.keys() else subst[tp] for tp in typelist]
        else:
            types = other_typelist
        for subsection in subsect_params:
            for parm_entry in subsection.entries_param:
                if parm_entry.match(types, int_type):
                    if from_wildtype is None:
                        if 'X' in parm_entry.types:
                            from_wildtype = True
                        else:
                            from_wildtype = False
                        new_params.append(parm_entry.params)
                        print(f'Fixing by analogy, using entry {str(parm_entry).strip()}')
                    elif (from_wildtype and 'X' in parm_entry.types) or (not from_wildtype
                                                                         and 'X' not in parm_entry.types):
                        new_params.append(parm_entry.params)
                        print(f'Fixing by analogy, using entry {str(parm_entry).strip()}')
                    if len(typelist) < 4:
                        break
            else:
                continue
            break
        return [y for x in new_params for y in x]
    
    def _add_ff_params_to_entry(self, entry, subsect_params):
        """
        Given a bonded term (e.g. "21     24     26    5") converts it to atomtypes,
        finds the respective FF parameters and adds them to the bonded entry
        :param entry: Entry, an EntryBonded instance to add FF params to
        :param subsect_params: list, SubsectionParam instances that hold all FF params
        :return: None
        """
        int_type = entry.interaction_type
        entry.read_types()
        for types, params, parmentry in zip([entry.types_state_a, entry.types_state_b],
                                            [entry.params_state_a, entry.params_state_b],
                                            [entry.params_state_a_entry, entry.params_state_b_entry]):
            wildcard_present = []
            non_wildcard_present = []
            for subsections in subsect_params:
                for parm_entry in subsections.entries_param:
                    if parm_entry.match(types, int_type):
                        is_wildcard = 'X' in parm_entry.types
                        if not wildcard_present and not is_wildcard:
                            params += parm_entry.params
                            parmentry.append(parm_entry)
                            non_wildcard_present += parm_entry.types
                        elif not wildcard_present and is_wildcard and not non_wildcard_present:
                            params += parm_entry.params
                            parmentry.append(parm_entry)
                            wildcard_present = parm_entry.types
                        elif wildcard_present and not is_wildcard:
                            raise RuntimeError("Wildcard ('X') entries were found prior to regular ones, please fix"
                                               "your FF parameters")
                        elif wildcard_present and is_wildcard:  # only add if multiple entries per given wildcard
                            if parm_entry.types == wildcard_present:
                                params += parm_entry.params
                                parmentry.append(parm_entry)
                            else:
                                pass
        if not entry.params_state_a and entry.subsection.header == 'dihedrals' \
                and any('DUM' in t for t in entry.types_state_a):
            entry.params_state_a = [0.0, 0.0, '1']
        if not entry.params_state_b and entry.subsection.header == 'dihedrals' and entry.types_state_b \
                and any('DUM' in t for t in entry.types_state_b):
            entry.params_state_b = [0.0, 0.0, '1']
        if entry.params_state_a and entry.subsection.header == 'dihedrals' and (not entry.params_state_b) \
                and entry.interaction_type in ('9', '4', '1'):
            if len(entry.params_state_a) > 3:
                assert len(entry.params_state_a) % 3 == 0
                leftover = entry.params_state_a[3:]
                entry.params_state_a = entry.params_state_a[:3]
                counter = 1
                while leftover:
                    new_entry = gml.EntryBonded(' '.join(str(x) for x in entry.content), self)
                    entry_location = entry.subsection.bkp_entries.index(entry)
                    entry.subsection.bkp_entries.insert(entry_location+counter, new_entry)
                    entry.subsection.bkp_entries[entry_location+counter].params_state_a = leftover[:3]
                    leftover = leftover[3:]
                    counter += 1
        if entry.params_state_a and entry.subsection.header == 'dihedrals' and entry.params_state_b \
                and entry.interaction_type in ('9', '4', '1'):
            if len(entry.params_state_a) > 3 or len(entry.params_state_b) > 3 \
                    or entry.params_state_a[2] != entry.params_state_b[2]:
                assert len(entry.params_state_a) % 3 == 0 and len(entry.params_state_b) % 3 == 0
                multiplicities_a = entry.params_state_a[2::3]
                multiplicities_b = entry.params_state_b[2::3]
                all_multiplicities = list(set(multiplicities_a + multiplicities_b))
                params_a = {entry.params_state_a[3*i+2]: entry.params_state_a[3*i:3*(i+1)]
                            for i in range(len(entry.params_state_a)//3)}
                params_b = {entry.params_state_b[3 * i + 2]: entry.params_state_b[3 * i:3 * (i + 1)]
                            for i in range(len(entry.params_state_b) // 3)}
                counter = 1
                m = all_multiplicities[-1]
                entry.params_state_a = params_a[m] if m in params_a.keys() else [0.0, 0.0, m]
                entry.params_state_b = params_b[m] if m in params_b.keys() else [0.0, 0.0, m]
                _ = all_multiplicities.pop()
                while all_multiplicities:
                    m = all_multiplicities[-1]
                    new_entry = gml.EntryBonded(' '.join(str(x) for x in entry.content), self)
                    entry_location = entry.subsection.bkp_entries.index(entry)
                    entry.subsection.bkp_entries.insert(entry_location + counter, new_entry)
                    entry.subsection.bkp_entries[entry_location + counter].params_state_a = params_a[m] \
                        if m in params_a.keys() else [0.0, 0.0, m]
                    entry.subsection.bkp_entries[entry_location + counter].params_state_b = params_b[m] \
                        if m in params_b.keys() else [0.0, 0.0, m]
                    counter += 1
                    _ = all_multiplicities.pop()
        if not entry.params_state_a and entry.subsection.header == 'dihedrals' and entry.params_state_b \
                and entry.interaction_type in ('9', '4', '1'):
            raise RuntimeError(f'Warning: in line {entry}, parameters were found for state A, but not for state B.'
                               f'Try to add parameters for types {entry.types_state_a}')
    
    def _check_parm_type(self):
        """
        Finds number codes for interaction type, e.g. CHARMM uses angletype '5' (urey-bradley)
        while Amber uses angletype '1' (simple harmonic); this will usually be one number
        but if subsections are merged, we can have mixed interaction types
        :return: list of strs, interaction types
        """
        types = []
        for entry in self.entries_bonded:
            if entry.interaction_type not in types:
                types.append(entry.interaction_type)
        return ['-1'] if not types else types

    def add_type_labels(self):
        for entry in self.entries_bonded:
            self._add_type_label(entry)

    @staticmethod
    def _add_type_label(entry):
        entry.read_types()
        entry.comment += " ; "
        entry.comment += " ".join(entry.types_state_a)
        if entry.types_state_b is not None:
            entry.comment += " ; "
            entry.comment += " ".join(entry.types_state_b)


class SubsectionParam(Subsection):
    """
    SubsectionParam contains force field parameters;
    should be included in SectionParam
    """
    n_atoms = {'pairtypes': 2, 'bondtypes': 2, 'constrainttypes': 2, 'angletypes': 3, 'dihedraltypes': 4,
               'nonbond_params': 2, 'defaults': 0, 'atomtypes': 1, 'implicit_genborn_params': 1, 'cmaptypes': 5}
    
    def __init__(self, content, section):
        super().__init__(content, section)
        self.atoms_per_entry = SubsectionParam.n_atoms[self.header]
        self.prmtypes = self._check_parm_type()
        self.label = '{}-{}'.format(self.header, self.prmtypes)
        self.ordering = {}
        if self.header == 'cmaptypes':
            self._process_cmap()
        
    def __repr__(self):
        if self.prmtypes[0] != '-1' or self.header not in ('atomtypes', 'implicit_genborn_params'):
            return "Subsection {} with interaction type {}".format(self.header, ' '.join(self.prmtypes))
        else:
            return "Subsection {}".format(self.header)
    
    def __add__(self, other):
        """
        Added for the purpose of merging subsections with
        identical headers
        :param other: other SubsectionParam instance
        :return: a new SubsectionParam instance resulting from the merger
        """
        if not isinstance(other, SubsectionParam):
            raise TypeError("{} is not a SubsectionParam instance".format(other))
        if self.header != other.header:
            raise TypeError("Cannot merge subsections with different headers: {} and {}".format(self.header,
                                                                                                other.header))
        return SubsectionParam(["[ {} ]\n".format(self.header)] + self.entries + other.entries, self.section)

    def get_entries_by_types(self, *types):
        found = []
        for entry in self.entries_param:
            if len(types) != len(entry.types):
                continue
            if all([query == tp for query, tp in zip(types, entry.types)]) \
                    or all([query == tp for query, tp in zip(types, entry.types[::-1])]):
                found.append(entry)
        if not found:
            raise RuntimeError(f"No entry found with types: {types}")
        return found

    @property
    def entries_param(self):
        return [e for e in self.entries if isinstance(e, gml.EntryParam)]

    def sort(self):
        """
        In case we want to sort entries after some are added at the end of the section
        :return: None
        """
        self.ordering = {tuple(str(entry).split()[:SubsectionParam.n_atoms[self.header]])
                         if isinstance(entry, gml.EntryParam) else str(entry): n
                         for n, entry in enumerate(self.entries)}
        self.entries.sort(key=self._sorting_fn)

    def find_used_ff_params(self):
        used_parm_entries = []
        used_atomtypes_a = {a.type for mol in self.section.top.molecules for a in mol.atoms}
        used_atomtypes_b = {a.type_b for mol in self.section.top.molecules for a in mol.atoms if a.type_b is not None}
        used_atomtypes = used_atomtypes_a.union(used_atomtypes_b)
        for entry in self.entries_param:
            if all(tp in used_atomtypes for tp in entry.types):
                used_parm_entries.append(entry.identifier)
        return used_parm_entries

    def _combine_entries(self, other_subsection: "gml.SubsectionParam", overwrite: bool = False):
        """
        When merging parameter sets, makes sure there are no duplications, and in these cases resolves conflicts
        :param other_subsection: the gml.Subsection instance from which entries will be added to merge files
        :param overwrite: bool, whether to accept new parameters in case of a mismatch
        :return: None
        """
        for entry in other_subsection.entries_param:
            existing = [e for e in self.entries_param if e.types == entry.types]
            if len(existing) == 1:
                ex_entry = existing[0]
                if ex_entry.params == entry.params and ex_entry.interaction_type == entry.interaction_type \
                        and ex_entry.modifiers == entry.modifiers:
                    continue
                else:
                    if overwrite:
                        ex_entry.params = entry.params
                        ex_entry.interaction_type = entry.interaction_type
                        ex_entry.modifiers = entry.modifiers
                    else:
                        print(f"Found mismatch between original entry {str(ex_entry)} and new {entry}, keeping the "
                              f"original; if this is not desired, set overwrite=True")
            elif len(existing) == 0:
                self.add_entry(entry)
            else:
                raise RuntimeError(f"Found multiple entries matching the types of {entry.types} in section"
                                   f"{self.header}, make sure this is intentional")

    def _sorting_fn(self, entry: "gml.Entry") -> int:
        """
        Comments should go first, then we sort based on first, second,
        ... column of the section
        :param entry: Entry, entry to be sorted
        :return: int, ordering number
        """
        types = tuple(str(entry).split()[:SubsectionParam.n_atoms[self.header]]) if isinstance(entry, gml.EntryParam) \
            else str(entry)
        line_number = self.ordering[types]
        return line_number if not hasattr(entry, 'types_state_a') or 'X' not in entry.types \
            else line_number + len(self.ordering.keys())
    
    def _check_parm_type(self):
        """
        Finds number code for interaction type, e.g. CHARMM uses angletype '5' (urey-bradley)
        while Amber uses angletype '1' (simple harmonic)
        :return: str, interaction type
        """
        if self.header not in SubsectionParam.n_atoms.keys() or self.header in ['atomtypes', 'implicit_genborn_params']:
            return '-1'
        npar = SubsectionParam.n_atoms[self.header]
        types = list({entry.content[npar] if len(entry.content) > npar else entry.interaction_type
                      for entry in self.entries_param})
        return '-1' if not types else types

    def get_opt_dih(self, types=False):
        dopts = [entry for entry in self.entries_param if 'DIHOPT' in entry.comment]
        if not types:
            return [e.params[x] for e in dopts for x in [0, 1]]
        else:
            return [(*e.types, e.params[2]) for e in dopts]

    def get_opt_dih_indices(self):
        dopts = [entry for entry in self.entries_param if 'DIHOPT' in entry.comment]
        self.section.top.add_ff_params()
        indices = set()
        for i in dopts:  # TODO so far only works for one (first) molecule
            for j in self.section.top.molecules[0].get_subsection('dihedrals').entries_bonded:
                if i in j.params_state_a_entry:
                    indices.add(j.atom_numbers)
        return list(indices)

    def set_opt_dih(self, values):
        dopts = [entry for entry in self.entries_param if 'DIHOPT' in entry.comment]
        for e, ang, k in zip(dopts, values[::2], values[1::2]):
            e.params[0] = ang
            e.params[1] = k

    def _remove_symm(self, prefix):
        for n in range(len(self.entries)):
            if n >= len(self.entries):
                break
            if isinstance(self.entries[n], gml.EntryParam) and any([x.startswith(prefix)
                                                                    for x in self.entries[n].types]):
                limit = min(8, len(self.entries)-n-1)
                next_other = [q for q in range(1, limit) if isinstance(self.entries[n+q], gml.EntryParam)
                              and self.entries[n].types != self.entries[n+q]]
                if not next_other:
                    continue
                else:
                    next_other = next_other[0]
                for other in range(n+next_other, len(self.entries)):
                    if other < len(self.entries) and isinstance(self.entries[other], gml.EntryParam) \
                            and any([x.startswith(prefix) for x in self.entries[other].types]):
                        if self.entries[n].types == self.entries[other].types[::-1] \
                                and self.entries[n].params == self.entries[other].params:
                            _ = self.entries.pop(other)
    
    def _process_cmap(self):
        """
        Reads a multiline entry from the [ cmaptypes ] section
        and converts it into an array that can be later properly printed
        :return: None
        """
        new_entries = []
        current = []
        for e in self.entries:
            if isinstance(e, gml.EntryParam):
                if e.content[-1].endswith('\\'):
                    current.extend([x.rstrip('\\') for x in e.content])
                else:
                    current.extend([x.rstrip('\\') for x in e.content])
                    new_entry = ' '.join(current)
                    new_entries.append(gml.EntryParam(new_entry, self, processed=True))
                    current = []
            else:
                new_entries.append(e)
        self.entries = new_entries
        
        
class SubsectionAtom(Subsection):
    """
    SubsectionAtom contains definitions of all atoms in the molecule;
    should be contained in SectionMol
    """
    def __init__(self, content, section):
        super().__init__(content, section)
        self.fstring = "{:>6}{:>11}{:>7}{:>7}{:>7}{:>7}{:>11}{:>11}   ; " + '\n'
        self.name_to_num, self.num_to_name, self.num_to_type, self.num_to_type_b = None, None, None, None

    @property
    def entries_atom(self):
        return [e for e in self.entries if isinstance(e, gml.EntryAtom)]

    def get_dicts(self, force_update=False):
        """
        dicts are not always needed and are costly to calculate,
        so only fill in the values when explicitly asked to
        :return: None
        """
        # TODO should they be @properties, called only when needed? EITHER WAY, REVISE
        if not self.name_to_num or force_update:
            self.name_to_num, self.num_to_name, self.num_to_type, self.num_to_type_b = self._mol_type_nums()

    def check_defined_types(self):
        defined_types = self.section.top.parameters.get_subsection('atomtypes')
        typelist = {ent.types[0] for ent in defined_types.entries_param}
        for atom in self.entries_atom:
            if atom.type not in typelist:
                print(f'Couldn\'t find definition of atom type {atom.type} (atom {atom.num} in molecule '
                      f'{self.section.mol_name}) in parameters')
            if atom.type_b:
                if atom.type_b not in typelist:
                    print(f'Couldn\'t find definition of atom type {atom.type_b} (atom {atom.num} in molecule '
                          f'{self.section.mol_name}) in parameters')

    def _mol_type_nums(self):
        """
        Provides bindings between atomnumber and atomtype
        and vice versa for each molecule identified in
        the topology
        :return: tuple of dicts, each dict contains molname:(type:num) and
        molname:(num:type) bindings
        """
        name_to_num, num_to_name, num_to_type, num_to_type_b = {}, {}, {}, {}
        for entry in self.entries_atom:
            name_to_num[(entry.resid, entry.atomname)] = entry.num
            num_to_name[entry.num] = entry.atomname
            num_to_type[entry.num] = entry.type
            num_to_type_b[entry.num] = entry.type_b if entry.type_b is not None else entry.type
        return name_to_num, num_to_name, num_to_type, num_to_type_b

    
class SubsectionHeader(Subsection):
    """
    SubsectionHeader contains the [ moleculetype ] section;
    should be contained in SectionMol
    """
    def __init__(self, content, section):
        super().__init__(content, section)
        self.molname = [a.content[0] for a in self.entries if a.content][0]

    @property
    def entries_bonded(self):
        return []

    @property
    def entries_atom(self):
        return []