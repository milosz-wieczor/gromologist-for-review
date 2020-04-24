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