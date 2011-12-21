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
