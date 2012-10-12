/*
 * COORD.cpp
 *
 *  Created on: Dec 21, 2011
 *      Author: asn
 */

#include "COORD.h"

COORD::COORD(PRMTOP* Mol, int as, int ae, int bs, int be, char* filename, bool gzipped) {
	this->astart = as;
	this->aend = ae;
	this->bstart = bs;
	this->bend = be;
	Energy = new ENERGY;
	if (gzipped){
		this->read_gzcrd(Mol, filename);
	}
	else {
		this->read_crd(Mol, filename);
	}
}

COORD::COORD(PRMTOP* Mol, int as, int ae, int bs, int be, char* filename) {
	this->astart = as;
	this->aend = ae;
	this->bstart = bs;
	this->bend = be;
	Energy = new ENERGY;
	this->read_netcdf(Mol, filename);
}


void COORD::read_crd(PRMTOP* Mol, char* filename) {
	step = 0;
	ifstream mdcrd(filename);
	getline(mdcrd, line);
	printf("# %s\n", line.c_str());
	printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");

	while (!mdcrd.eof()){
		current_crd.clear();
		for (int i=0; i < Mol->N; i++) {
			xyz.clear();
			for (int j=0; j < 3; j++) {
				mdcrd >> pos;
				xyz.push_back(pos);
			}
			current_crd.push_back(xyz);
		}
		if (Mol->ifbox > 0){
			for (int i=0; i<3; i++){
				mdcrd >> pos;	// BOX INFO
			}
		}
		step++;
		printf(" %12d ", step);
		Energy->compute_nb2(Mol, current_crd, this->astart, this->aend, this->bstart, this->bend);
	}
}

void COORD::read_gzcrd(PRMTOP* Mol, char* filename) {
	step = 0;
	igzstream mdcrd(filename);
	getline(mdcrd, line);
	printf("# %s\n", line.c_str());
	printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");

	while (!mdcrd.eof()){
		current_crd.clear();
		for (int i=0; i < Mol->N; i++) {
			xyz.clear();
			for (int j=0; j < 3; j++) {
				mdcrd >> pos;
				xyz.push_back(pos);
			}
			current_crd.push_back(xyz);
		}
		if (Mol->ifbox > 0 ){
			for (int i=0; i< 3; i++){
				mdcrd >> pos; 	// Box Info
			}
		}
		step++;
		printf(" %12d ", step);
		Energy->compute_nb2(Mol, current_crd, this->astart, this->aend, this->bstart, this->bend);
	}
}

COORD::~COORD() {
	current_crd.clear();
	xyz.clear();
}

void COORD::read_netcdf(PRMTOP* Mol, char* filename){
	NcFile nc_mdcrd(filename, NcFile::ReadOnly);

	if (!nc_mdcrd.is_valid()){
		printf("$ Could not open trajectory file %s. Please check.\n", filename);
	}

	NcDim* FrameDim = nc_mdcrd.get_dim("frame");
	int size = FrameDim->size();
	printf("# NetCDF Frame dimension: %d\n", size);

	NcDim* NDim = nc_mdcrd.get_dim("atom");
	const int Ndim = NDim->size();
	if (Ndim != Mol->N){
		printf("# Mismatch among number of atoms in PRMTOP (%d) and NETCDF (%d) files. Please check.\n", Mol->N, Ndim);
		exit(1);
	}
	else {
		printf("# NetCDF number of atoms: %d\n", Ndim);
	}

	NcVar* nc_Coordinates = nc_mdcrd.get_var("coordinates");
	double coords[Ndim][3];

	printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");

	for (int frame=1; frame <= size; frame++){
		nc_Coordinates->get(&coords[0][0], 1, Ndim, 3);
// test
		current_crd.clear();
		for (int i=0; i< Ndim; i++){
			for (int j=0; j< 3; j++){
				xyz.push_back(coords[i][j]);
			}
			printf("%d --> %7.3f  %7.3f  %7.3f\n", i, xyz[0], xyz[1], xyz[2]);
			current_crd.push_back(xyz);
			xyz.clear();
		}
		printf(" %12d ", frame);
//		Energy->compute_nb2(Mol, coords, this->astart, this->aend, this->bstart, this->bend);
		Energy->compute_nb2(Mol, current_crd, this->astart, this->aend, this->bstart, this->bend);
		nc_Coordinates->set_cur(frame, 0, 0);
	}
}
