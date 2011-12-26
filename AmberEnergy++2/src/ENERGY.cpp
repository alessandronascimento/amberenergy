/*
 * ENERGY.cpp
 *
 *  Created on: Dec 21, 2011
 *      Author: asn
 */

#include "ENERGY.h"

ENERGY::ENERGY() {
}

double ENERGY::compute_r(vector<double> i, vector<double> j) {
	x= j[0]-i[0];
	y= j[1]-i[1];
	z= j[2]-i[2];
	r= sqrt((x*x)+(y*y)+(z*z)) ;
	return (r) ;
}

double ENERGY::compute_nb2 (PRMTOP* Mol, vector<vector<double> > coords, int astart, int aend, int bstart, int bend) {

	elec=0.00;
	vdw=0.00;
	for (int i=astart; i<=aend; i++) {
		this->get_excluded_atoms(Mol, i);

		for (int j=bstart; j<=bend; j++) {

			this->excluded=false;
			for(unsigned a=0; a< this->excluded_atoms.size(); a++){
				if (excluded_atoms[a]-1 == j){
					excluded=true;
				}
			}

			if ((j < astart or j > aend) and (not this->excluded)){

				r = this->compute_r(coords[i], coords[j]);
				r2 = r*r;
				elec += ((Mol->charges[i]*Mol->charges[j])/(r));
				iaci = Mol->Natomtypes*(Mol->atom_types_index[i]-1);
				ic = Mol->ico[(iaci + Mol->atom_types_index[j])-1];

#ifdef DEBUG
				printf("ico for atoms %s and %s: %d\n", Mol->atomnames[i].c_str(), Mol->atomnames[j].c_str(), ic);
#endif

				if (ic > 0) { 		// Use 12 - 6 LJ Potencial
					vdw += (((Mol->CN1[ic-1])/(r2*r2*r2*r2*r2*r2))-((Mol->CN2[ic-1])/(r2*r2*r2)));
#ifdef DEBUG
					printf("CN1: %10.4f    CN2: %10.4f\n", Mol->CN1[ic-1], Mol->CN2[ic-1]);
#endif
				}

				else {			// Use 10 - 12 LJ Potencial
					vdw += (((Mol->ASOL[-(ic-1)])/(r2*r2*r2*r2*r2*r2)) - ((Mol->BSOL[-(ic-1)])/(r2*r2*r2*r2*r2)));
#ifdef DEBUG
					printf("ASOL: %10.4f    BSOL: %10.4f\n", Mol->ASOL[-(ic-1)], Mol->CN2[-(ic-1)]);
#endif
				}
			}
		}
	}
	printf("%12.4f %12.4f %12.4f\n", elec, vdw, (elec+vdw));
	return (elec+vdw);
}

void ENERGY::get_excluded_atoms(PRMTOP* Mol, int atom) {
	exc_atoms=0;
	exc_start=0;
	exc_end=0;

	exc_atoms=Mol->number_excluded_atoms[atom];

	if (exc_atoms > 0){
		exc_start=0;
		for (int a=0; a<atom; a++){
			exc_start += Mol->number_excluded_atoms[a];
		}

		exc_end=exc_start+exc_atoms-1;
	}

	this->excluded_atoms.clear();

	for (int a=exc_start; a<=exc_end; a++){
		excluded_atoms.push_back(Mol->excluded_atoms_list[a]);
	}

#ifdef DEBUG
	printf ("Atom %d has    %d excluded atoms (%d to %d).\n", i, excluded_atoms.size(), exc_start, exc_end);
	for (int a=exc_start; a<=exc_end; a++){
		printf(" %d", Mol->excluded_atoms_list[a]-1);
	}
	printf("\n");
#endif

}


ENERGY::~ENERGY() {
}
