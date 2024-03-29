/*
 * COORD.h
 *
 *  Created on: Dec 21, 2011
 *      Author: asn
 */

#ifndef COORD_H_
#define COORD_H_

#include "PRMTOP.h"
#include "ENERGY.h"
#include "gzstream.h"
#include <fstream>
#include <stdio.h>
#include <string>
//#include <netcdf.h>
#include <netcdfcpp.h>

class COORD {
public:
	COORD(PRMTOP* Mol, int as, int ae, int bs, int be, char* filename, bool gzipped);
	COORD(PRMTOP* Mol, int as, int ae, int bs, int be, char* filename, int mode);
    COORD(PRMTOP* Mol, vector<int> residues_atoms, vector<int> ligand_atoms, char* filename, int mode);
	int astart, aend, bstart, bend;
    vector<int> receptor_atom_chuncks;
    vector<int> ligand_atom_chuncks;
	vector<vector<double> > current_crd;
	vector<double> xyz;
	string line;
	double pos;
	unsigned step;
	ENERGY* Energy;
	virtual ~COORD();
	void read_crd(PRMTOP* Mol, char* filename);
	void read_gzcrd(PRMTOP* Mol, char* filename);
	void read_netcdf(PRMTOP* Mol, char* filename);
    void read_netcdf(PRMTOP* Mol, char* filename, vector<int> rec, vector<int> lig);
	void read_dcd(PRMTOP* Mol, char* filename);
	void parse_binary_int(int* i, FILE* file);
};

#endif /* COORD_H_ */
