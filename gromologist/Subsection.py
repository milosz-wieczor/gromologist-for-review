import gromologist as gml
from copy import deepcopy


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
        self.header = content[0].strip().strip('[]').strip()
        if ';' in self.header:
            pos = self.header.index(';')
            self.header = self.header[:pos].strip().strip('[]').strip()
        self.write_header = self.header if self.header != 'impropers' else 'dihedrals'
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
        :param new_entry: str, entry to be added
        :param position: where to add the entry (None is at the end)
        :return: None
        """
        if position:
            position = int(position)
            self.entries.insert(position, new_entry)
        else:
            self.entries.append(new_entry)
    
    def add_entries(self, new_entries_list, position=None):
        """
        Adds multiple entries to the subsection, either at the end
        or in a specified position
        :param new_entries_list: list of str, entries to be added
        :param position: where to add the entries (None is at the end)
        :return: None
        """
        if position:
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
        
        
class SubsectionBonded(Subsection):
    """
    SubsectionBonded contains a subsection with entries corresponding to bonded terms,
    e.g., bonds or dihedrals; should be included in SectionMol
    """
    n_atoms = {'bonds': 2, 'pairs': 2, 'angles': 3, 'dihedrals': 4, 'impropers': 4,
               'cmap': 5, 'settles': 2, 'exclusions': 2, 'position_restraints': 1}
    
    def __init__(self, content, section):
        super().__init__(content, section)
        self.atoms_per_entry = SubsectionBonded.n_atoms[self.header]
        self.prmtype = self._check_parm_type()
        self.label = '{}-{}'.format(self.header, self.prmtype)
        self.fstring = "{:5} " * (SubsectionBonded.n_atoms[self.header] + 1) + '\n'
    
    def __repr__(self):
        return "Subsection {} with interaction type {}".format(self.header, self.prmtype)
    
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
    
    def add_ff_params(self):
        matchings = {'bonds': 'bondtypes', 'angles': 'angletypes', 'dihedrals': 'dihedraltypes',
                     'impropers': 'dihedraltypes'}
        sections = [sect for sect in self.section.top.sections if sect.name == 'Parameters']
        subsect_params = [sub for sect in sections for sub in sect.subsections
                          if sub.header == matchings[self.header]]
        self.bkp_entries = self.entries[:]
        for entry in self.entries:
            if isinstance(entry, gml.EntryBonded):
                self._add_ff_params_to_entry(entry, subsect_params)
        self.entries = self.bkp_entries[:]
    
    def _add_ff_params_to_entry(self, entry, subsect_params):
        """
        Given a bonded term (e.g. "21     24     26    5") converts it to atomtypes,
        finds the respective FF parameters and adds them to the bonded entry
        :param entry: Entry, an EntryBonded instance to add FF params to
        :param subsect_params: list, SubsectionParam instances that hold all FF params
        :return: None
        """  # TODO add some info to comments
        int_type = entry.interaction_type
        entry.read_types()
        wildcard_present = []
        non_wildcard_present = []
        for subsections in subsect_params:
            for parm_entry in [e for e in subsections if isinstance(e, gml.EntryParam)]:
                if parm_entry.match(entry.types_state_a, int_type):
                    is_wildcard = 'X' in parm_entry.types
                    if not wildcard_present and not is_wildcard:
                        entry.params_state_a += parm_entry.params
                        non_wildcard_present += parm_entry.types
                    elif not wildcard_present and is_wildcard and not non_wildcard_present:
                        entry.params_state_a += parm_entry.params
                        wildcard_present = parm_entry.types
                    elif wildcard_present and not is_wildcard:
                        raise RuntimeError("Wildcard ('X') entries were found prior to regular ones, please fix"
                                           "your FF parameters")
                    elif wildcard_present and is_wildcard:  # only add if there are multiple entries per given wildcard
                        if parm_entry.types == wildcard_present:
                            entry.params_state_a += parm_entry.params
                        else:
                            pass
        if entry.params_state_a and entry.subsection.header == 'dihedrals':
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
    
    def _check_parm_type(self):
        """
        Finds number code for interaction type, e.g. CHARMM uses angletype '5' (urey-bradley)
        while Amber uses angletype '1' (simple harmonic)
        :return: str, interaction type
        """
        for entry in self:
            if isinstance(entry, gml.EntryBonded):
                return entry.interaction_type
        return '0'


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
        self.prmtype = self._check_parm_type()
        self.label = '{}-{}'.format(self.header, self.prmtype)
        if self.header == 'cmaptypes':
            self._process_cmap()
        
    def __repr__(self):
        if self.prmtype != '0' or self.header not in ('atomtypes', 'implicit_genborn_params'):
            return "Subsection {} with interaction type {}".format(self.header, self.prmtype)
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
    
    def _check_parm_type(self):
        """
        Finds number code for interaction type, e.g. CHARMM uses angletype '5' (urey-bradley)
        while Amber uses angletype '1' (simple harmonic)
        :return: str, interaction type
        """
        if self.header not in SubsectionParam.n_atoms.keys():
            return '0'
        npar = SubsectionParam.n_atoms[self.header]
        for entry in self:
            if len(entry.content) > npar and isinstance(entry, gml.EntryParam):
                return entry.content[npar]
        return '0'
    
    def _process_cmap(self):
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
        self.nat, self.charge = None, None
        self.calc_properties()
        self.name_to_num, self.num_to_name, self.num_to_type, self.num_to_type_b = None, None, None, None
    
    def calc_properties(self):
        self.nat = self.section.natoms = self._calc_nat()
        self.charge = self.section.charge = self._calc_charge()
        
    def _calc_charge(self):
        """
        Calculates total charge of the molecule
        :return: float, total charge
        """
        total_charge = 0
        for entry in self.entries:
            if isinstance(entry, gml.EntryAtom):
                total_charge += entry.charge
        return total_charge
    
    def _calc_nat(self):
        return len([e for e in self.entries if isinstance(e, gml.EntryAtom)])

    def _get_dicts(self):
        """
        dicts are not always needed and are costly to calculate,
        so only fill in the values when explicitly asked to
        :return: None
        """
        if not self.name_to_num:
            self.name_to_num, self.num_to_name, self.num_to_type, self.num_to_type_b = self._mol_type_nums()

    def _mol_type_nums(self):
        """
        Provides bindings between atomnumber and atomtype
        and vice versa for each molecule identified in
        the topology
        :return: tuple of dicts, each dict contains molname:(type:num) and
        molname:(num:type) bindings
        """
        name_to_num, num_to_name, num_to_type, num_to_type_b = {}, {}, {}, {}
        for entry in self:
            if isinstance(entry, gml.EntryAtom):
                name_to_num[entry.atomname] = entry.num
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
