# iPot

**iPot** (from **i**nverse **Pot**ts) - a set of programs to calculate inter- and intra- protein contact energies using 
atom-atom and residue-residue statistical contact potentials derived from the Potts model [1].

```
aace18    - atomic contact energies (18 atom types from Ref. [2])
aace20    - atomic contact energies (20 types, all heavy atoms in a residue belong to one type)
aace167   - atomic contact energies (167 types, all heavy atoms in a residue belong to different types)
rrce20    - contact energies between residue centroids
```

`aace18` with default parameters should be sufficient for most applications.

## Installation

#### Program download and compilation

```
git clone https://github.com/gjoni/aace-rrce
cd aace-rrce
make
```

#### Download contact potentials database (optional)
```
wget http://vakser.bioinformatics.ku.edu/resources/ipot/download/ipot_data.tar.gz
tar xf ipot_data.tar.gz
```

## Usage

By the example of `aace18` (other programs follow similar logic):

```
Usage:   ./aace18 [-option] [argument]

Options:  -r receptor.pdb \  # PDB file with receptor's coordinates
          -l ligand.pdb \    # (optional) PDB file with ligand's coordinates
          -t table.txt \     # (optional) aace18 contact potential table 
          -d dmax            # (optional) contact distance, default dmax=6.8A
```

##### Intrachain energy:
```
./aace18 -r example/1i2m_u1.pdb
```

##### Interchain energy:
```
./aace18 -r example/1i2m_u1.pdb -l example/1i2m_u2.pdb
```

##### Interchain energy using non-default energy table:
```
./aace18 -r example/1i2m_u1.pdb -l example/1i2m_u2.pdb -t /path/to/database/ipot_data/AACE18/table.8.0A_k4 -d 8.0
```


## Acknowledgements

 - [kdtree library](https://github.com/jtsiomb/kdtree) by John Tsiombikas

## Links

 - [Vakser Lab](http://vakser.compbio.ku.edu/main/)

## References
[1] I Anishchenko, PJ Kundrotas, IA Vakser. Contact potential for structure prediction 
of proteins and protein complexes from Potts model. (2018) In preparation

[2] C Zhang, G Vasmatzis, JL Cornette, C DeLisi. Determination of atomic 
desolvation energies from the structures of crystallized proteins. (1997) 
J Mol Biol. 267(3): 707-6
