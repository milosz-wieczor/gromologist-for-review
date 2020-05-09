# Gromologist

Gromologist is a package designed to facilitate handling, editing and manipulating GROMACS topology files (.top and .itp).

## Installation

After cloning to a local repository (`git clone https://gitlab.com/KomBioMol/gromologist.git`),
the package can be installed into Python by typing `pip install .` in the main Gromologist directory.
If you're using Anaconda, the same will work with `path/to/anaconda/bin/pip`.

## Usage

`Top` and `Pdb` are the core classes of the library, and are supposed to provide representation
for topology and structure objects, respectively. To initialize them, a path to the file
should be passed to the constructor:

```
>>> from gromologist import Top, Pdb
>>> t = Top('md/topol.top')
>>> p = Pdb('md/conf.pdb')
```

Since all .itp files are by default included into the `Top` object, sometimes it is
necessary to specify a custom path to Gromacs database:

```
>>> t = Top('md/topol.top', gmx_dir='/home/user/gromacs/share/gromacs/top')
```

Alternatively, `Top` can be initialized with both paths, or `Pdb` can be supplied later.
Note that one `Top` instance can only hold one `Pdb` instance at a time.

```
>>> t = Top('md/topol.top', pdb='md/conf.pdb')
>>> t.pdb
PDB file md/conf.pdb with 100 atoms
>>> t.add_pdb('md/other.pdb')
>>> t.pdb
PDB file md/other.pdb with 105 atoms
```

If `Pdb` is bound to `Top`, a number of diagnostic and fixing options are available,
including name consistency checks:

```
>>> t.check_pdb()
Check passed, all names match
```

as well as renumbering, automatic addition of chains and fixing ordering issues in PDB
using topology ordering as a template:

```
>>> t.pdb.renumber_all()
>>> t.pdb.add_chains()
>>> t.pdb.match_order_by_top_names()
```

With `Pdb.select_atoms()`, selections can be made in a VMD-like manner:

```
>>> t.pdb.select_atoms('name CA and (resname ASP or chain B)')
[5, 60, 72, 88]
```

After changes have been made, modified files can be saved:

```
>>> t.save_top('md/new_topol.top')
>>> t.pdb.save_pdb('md/new_conf.pdb')
```

If the topology contains too many parts irrelevant to the system at hand,
a leaner version can be produced that lacks unused molecule definitions:

```
>>> t.clean_sections()
```

To save a 'lightweight' .top file with all contents split into separate .itp files, 
use the `split` parameter of `Top.save_top`:

```
>>> t.save_top('md/new_topol.top', split=True) 
```

### Editing topologies

Let's start with a generic topology file:

```
>>> t = Top('md/topol.top')
```

##### Adding bonds within or between molecules

One useful application of Gromologist is adding bonds (and, automatically, other bonded terms)
either within a molecule or between them:

```
>>> protein = t.get_molecule("Protein")
>>> ligand = t.get_molecule("Other")
>>> t.list_molecules()
Protein                      1
Other                        1
>>> protein.merge_two(ligand, anhor_own=5, anchor_other=1)
>>> t.list_molecules()
Protein                      1
```

The above script merges Protein and Other into a single Protein molecule, adding a bond
between atom 5 of Protein and atom 1 of Other (here, indices are 1-based, corresponding
to numbering in .itp files).

To add a bond within a single e.g. Protein molecule, one can use `protein.merge_two(protein, anhor_own=2, anchor_other=3)`
or, more simply, `protein.add_bond(5,3)`.

##### Adding and removing atoms while maintaining ordered numbering

When an atom is removed, other atom numbers are modified accordingly, something that has to be
considered when removing multiple atoms. For instance, one can remove the first three atoms
in the following manner:

```
>>> protein.del_atom(1)
>>> protein.del_atom(1)
>>> protein.del_atom(1)
```

Note that all parameters involving this atom are automatically removed as well.

To add an atom, one should specify its desired placement within the molecule, and at least 
a minimal set of parameters:

```
>>> protein.add_atom(atom_number=20, atom_name="CA", atom_type="CT")
By default, atom will be assigned to residue MET1. Proceed? [y/n]
y
>>> protein.add_bond(20,8)
```

If residue data is not specified, Gromologist will attempt to guess the residue based on
neighboring atoms.

##### Adding alchemical B-states

To generate alchemical states for a subset of atoms, one can use `gen_state_b`:

```
>>> protein.gen_state_b(atomtype='CT',new_type="CE")
```

The arguments for `gen_state_b` are divided into two subgroups:

 + `atomname`, `resname`, `resid` and `atomtype` behave as selectors, allowing to specify
 one or many atoms that should have its B-type specified;
 + `new_type`, `new_charge` and `new_mass` act as setters, allowing to specify the values
 in the B-state.
 
If the selectors point to multiple atoms (e.g. `atomtype=CT` selects all atoms with type CT),
all will be modified as specified. In turn, if a given setter is not specified, the respective 
value will be identical to that for state A.

##### Duplicating types

Often it's useful to duplicate an atomtype exactly, i.e., assign it a different name while
retaining all bonded and nonbonded parameters of the original. This can be done easily with:

```
>>> params = t.parameters
>>> params.clone_type("CT", prefix="Y")
```

This will create a type "YCT" that shares all properties with "CT" but can be modified independently.

##### Explicitly listing parameters in topology 

To explicitly include all parameters in sections `[ bonds ]`, `[ angles ]` and `[ dihedrals ]`,
one can use:

```
>>> t.add_ff_params()
>>> t.save_top('with_params.top')
```

### Editing structures

Let's start by reading a PDB file:

```
>>> p = PDB('md/other.pdb')
```

##### Adding atoms along a vector specified by other atoms

To add e.g. a hydrogen atom to an existing atom CB, with a bond length of 1 A in the direction
specified by a vector from atom C to atom CA, one can use:

```
>>> p.insert_atom(len(p.atoms), base_atom=p.atoms[11], atomname='HX')
>>> p.reposition_atom_from_hook('name HX', 'name CB', 1.0, p1_sel='name C', p2_sel='name CA')
```

All the selections should be unique (corresponding to a single atom), and combinations can be
used like in VMD, e.g. `name CB and resid 10 and chain A`.

##### Filling beta-values with custom data (for visualization)

To use the PDB's beta column for color-mapping of observables e.g. in VMD, use the following:

```
>>> ca_atom_indices = p.select_atoms('name CA')
>>> p.fill_beta(per_residue_data, serials=[x+1 for x in ca_atom_indices])
>>> p.save_pdb('with_betas.pdb')
```

By adding the `smooth=...` parameter to `Pdb.fill_beta`, data can be spatially smoothed
using a Gaussian kernel with a specified standard deviation (in A).

##### Creating new PDB as a subset of existing one

To choose and save e.g. only the DNA atoms from a protein-DNA complex, use:

```
>>> dna =  Pdb.from_selection(p, 'dna')
>>> dna.save_pdb('dna.pdb')
```
