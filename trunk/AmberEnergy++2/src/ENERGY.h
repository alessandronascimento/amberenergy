/*
 * ENERGY.h
 *
 *  Created on: Dec 21, 2011
 *      Author: asn
 */

#ifndef ENERGY_H_
#define ENERGY_H_

#include "PRMTOP.h"
#include <math.h>
#include <stdio.h>

class ENERGY {
public:
	double elec;
	double vdw;
	double x, y, z;
	double r, r2;
	int iaci, ic;
	ENERGY();
	virtual ~ENERGY();

	double compute_r(vector<double> i, vector<double> j);
	double compute_nb2 (PRMTOP* Mol, vector<vector<double> > coords, int astart, int aend, int bstart, int bend);
};

#endif /* ENERGY_H_ */
