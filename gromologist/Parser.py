from .Entries import EntryAtom


class SelectionParser:
    def __init__(self, master):
        self.master = master
        try:
            self.nat = self.master.natoms
        except AttributeError:
            raise TypeError("Can only parse selections with PDB or topology data")
        
    @staticmethod
    def ispdb(obj) -> bool:
        try:
            _ = obj._cryst_format
        except AttributeError:
            return False
        else:
            return True
    
    def __call__(self, selection_string: str) -> list:
        while "  " in selection_string:
            selection_string = selection_string.replace("  ", " ")
        selection_string = selection_string.strip()
        protein_selection = "resname ALA ACE CYS ASP ASPP GLU GLUP PHE GLY HIS HID HIE HSD HSE ILE LYS LEU MET " \
                            "NME NMA ASN PRO GLN ARG SER THR VAL TYR TRP"
        dna_selection = "resname DA DG DC DT DA5 DG5 DC5 DT5 DA3 DG3 DC3 DT3"
        rna_selection = "resname RA RG RC RU RA5 RG5 RC5 RU5 RA3 RG3 RC3 RU3 A G C U A3 A5 G3 G5 C3 C5 U 3 U5"
        solvent_selection = "resname HOH TIP3 SOL OPC SPC WAT K CL NA POT K+ NA+ CLA CL-"
        all_selection = "serial < 10000000"
        noh_selection = "not element H"
        selection_string = selection_string.replace('solvent', solvent_selection)
        selection_string = selection_string.replace('water', 'resname HOH TIP3 SOL OPC SPC WAT')
        selection_string = selection_string.replace('protein', protein_selection)
        selection_string = selection_string.replace('nucleic', '(dna or rna)')
        selection_string = selection_string.replace('noh', noh_selection)
        selection_string = selection_string.replace('dna', dna_selection)
        selection_string = selection_string.replace('rna', rna_selection)
        selection_string = selection_string.replace('all', all_selection)
        return sorted(list(self._select_set_atoms(selection_string)))
    
    def _select_set_atoms(self, selection_string: str) -> set:
        """
        Main recursive fn taking care of stratifying the input
        :param selection_string: str, selection string
        :return: set of int, output atom indices
        """
        assert isinstance(selection_string, str)
        selection_string = selection_string.strip()
        parenth_ranges, operators = self._parse_sel_string(selection_string)
        # if there's just one parenthesis around the full phrase, let's remove it:
        while len(parenth_ranges) == 1 and selection_string[0] == '(' and selection_string[-1] == ')':
            selection_string = selection_string[1:-1]
            parenth_ranges, operators = self._parse_sel_string(selection_string)
        # simple phrase with no parenthesis and no operators, except for "not":
        if not parenth_ranges and not operators:
            if selection_string.strip().startswith("not "):
                return self._find_atoms(selection_string.lstrip()[4:], rev=True)
            else:
                return self._find_atoms(selection_string)
        # parenthesis but no operators, except for "not":
        elif parenth_ranges and not operators:
            if selection_string.strip().startswith("not "):
                set_all = {n for n in range(self.nat)}
                return set_all.difference(self._select_set_atoms(selection_string.strip()[4:].strip('()')))
            else:
                return self._select_set_atoms(selection_string.strip().strip('()'))
        else:  # we have operators and (possibly) parentheses
            if not parenth_ranges:
                first_op_borders = operators[0]
                first_op = selection_string[operators[0][0]:operators[0][1]].strip()
            else:
                first_op_borders = [opb for opb in operators if all([opb[0] not in range(*par) for par in parenth_ranges])][0]
                first_op = selection_string[first_op_borders[0]:first_op_borders[1]].strip()
            if first_op == "and":
                return self._select_set_atoms(selection_string[:first_op_borders[0]]) \
                    .intersection(self._select_set_atoms(selection_string[first_op_borders[1]:]))
            elif first_op == "or":
                return self._select_set_atoms(selection_string[:first_op_borders[0]]) \
                    .union(self._select_set_atoms(selection_string[first_op_borders[1]:]))
            elif first_op.startswith("same"):
                return self.master.same_residue_as(self._select_set_atoms(selection_string[first_op_borders[1]:]))
            elif first_op.startswith("within") or first_op.startswith("pbwithin"):
                if not self.ispdb(self.master):
                    raise ValueError("the within keyword only works for structural data, not topology")
                nopbc = False if first_op == "pbwithin" else True
                within = float([x for x in first_op.split() if x.isnumeric()][0])
                return self.master.within(self._select_set_atoms(selection_string[first_op_borders[1]:]), within, nopbc=nopbc)
    
    @staticmethod
    def _parse_sel_string(selection_string: str) -> (list, list):
        parenth_ranges = []
        operators = []
        opened_parenth = 0
        beginning = -1
        for nc, char in enumerate(selection_string):
            if char == '(':
                opened_parenth += 1
                if beginning == -1:
                    beginning = nc
            elif char == ')':
                opened_parenth -= 1
                end = nc
                if opened_parenth == 0:
                    parenth_ranges.append((beginning, end))
                    beginning = -1
            if opened_parenth < 0:
                raise ValueError("Improper use of parentheses in selection string {}".format(selection_string))
            if selection_string[nc:nc + 5] == " and " and opened_parenth == 0:
                operators.append((nc, nc + 5))
            if selection_string[nc:nc + 4] == " or " and opened_parenth == 0:
                operators.append((nc, nc + 4))
            if (selection_string[nc:nc + 7] == "within " or selection_string[nc:nc + 9] == "pbwithin ") and opened_parenth == 0:
                ending = selection_string.find(" of")
                operators.append((nc, ending + 4))
            if selection_string[nc:nc+16] == "same residue as ":
                operators.append((nc, nc + 16))
        if opened_parenth != 0:
            raise ValueError("Improper use of parentheses in selection string {}".format(selection_string))
        return parenth_ranges, operators
    
    def _find_atoms(self, sel_string: str, rev: bool = False) -> set:
        chosen = []
        keyw = sel_string.split()[0]
        if self.ispdb(self.master):
            matchings = {"name": "atomname", "resid": "resnum", "resnum": "resnum", "element": "element",
                         "chain": "chain", "resname": "resname", "serial": "serial"}
        else:
            matchings = {"name": "atomname", "resid": "resid", "resnum": "resid", "mass": "mass", "element": "element",
                         "resname": "resname", "serial": "num", "type": "type", "molecule": "molname"}
        try:
            vals = {int(x) for x in sel_string.split()[1:]}
        except ValueError:
            if " to " in sel_string:
                beg = int(sel_string.split()[1])
                end = int(sel_string.split()[3])
                vals = range(beg, end + 1)
            elif "<=" in sel_string:
                beg = 0
                end = int(sel_string.split('<=')[1].strip())
                vals = range(beg, end + 1)
            elif ">=" in sel_string:
                beg = int(sel_string.split('>=')[1].strip())
                end = 10**7
                vals = range(beg, end)
            elif "<" in sel_string:
                beg = 0
                end = int(sel_string.split('<')[1].strip())
                vals = range(beg, end)
            elif ">" in sel_string:
                beg = int(sel_string.split('>')[1].strip())
                end = 10**7
                vals = range(beg+1, end)
            else:
                vals = set(sel_string.split()[1:])
        atomlist = self.master.atoms
        for n, a in enumerate(atomlist):
            if not rev:
                if a.__getattribute__(matchings[keyw]) in vals:
                    chosen.append(n)
            else:
                if a.__getattribute__(matchings[keyw]) not in vals:
                    chosen.append(n)
        return set(chosen)
