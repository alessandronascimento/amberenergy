/*
 * ENERGY.cpp
 *
 *  Created on: Dec 21, 2011
 *      Author: asn
 */

#include "ENERGY.h"

ENERGY::ENERGY() {
	// TODO Auto-generated constructor stub

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
		for (int j=bstart; j<=bend; j++) {
			r = this->compute_r(coords[i], coords[j]);
			r2 = r*r;
			elec += ((Mol->charges[i]*Mol->charges[j])/(r));

			iaci = Mol->Natomtypes*(Mol->atom_types_index[i]-1);
			ic = Mol->ico[(iaci+Mol->atom_types_index[j])-1];

#ifdef DEBUG
			printf("ico for atoms %s and %s: %d        ", Mol->atomnames[i].c_str(), Mol->atomnames[j].c_str(), ic);
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
	printf("%12.4f %12.4f %12.4f\n", elec, vdw, (elec+vdw));
	return (elec+vdw);
}


ENERGY::~ENERGY() {
	// TODO Auto-generated destructor stub
}
