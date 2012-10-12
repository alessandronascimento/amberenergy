/*
 * ENERGY.h
 *
 *  Created on: Dec 21, 2011
 *      Author: asn
 */

#ifndef ENERGY_H_
#define ENERGY_H_

#include "PRMTOP.h"
#include <vector>
#include <math.h>
#include <stdio.h>

class ENERGY {
public:
	double elec;
	double vdw;
	double x, y, z;
	double r, r2;
	int iaci, ic;
	int exc_atoms;
	int exc_start;
	int exc_end;
	vector<int> excluded_atoms;
	bool excluded;
	ENERGY();
	virtual ~ENERGY();

	double compute_r(vector<double> i, vector<double> j);
	double compute_r(double i[], double j[]);
	double compute_nb2 (PRMTOP* Mol, vector<vector<double> > coords, int astart, int aend, int bstart, int bend);
	double compute_nb2 (PRMTOP* Mol, double coords[][3], int astart, int aend, int bstart, int bend);
	void get_excluded_atoms(PRMTOP* Mol, int atom);
};

#endif /* ENERGY_H_ */
