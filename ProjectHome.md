This is a really simple code to:
  * read a file containing the topologies/parameters of a system (prmtop);
  * read a MM simulation trajectory file in amber format;
  * compute interaction energies in the trajectory.

It was originally developed to compute interaction energies between a ligand and a given residue in the active site, for example. However, but might be useful to evaluate interaction energies between two residues in a protein or even between one residue/ligand and the whole environment (i.e., all other atoms in the system).

In order to compile the provided source code, the user will need:
  * a c++ compiler (works fine with gnu g++, intel, MinGW);
  * the the gzstream library (which requires the zlib library). This libray adds the support to read the gzipped trajectory files. This library can be obtained from http://www.cs.unc.edu/Research/compgeom/gzstream.
  * The netcdf library ( http://www.unidata.ucar.edu/software/netcdf );


Supported trajectory formats (so far):
  * Amber coordinates (mdcrd) files;
  * NetCDF trajactory format;
  * CHARMM/NAMD/X-PLOR DCD file format.