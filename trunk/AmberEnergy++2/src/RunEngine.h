/*
 * RunEngine.h
 *
 *  Created on: Dec 21, 2011
 *      Author: asn
 */

#ifndef RUNENGINE_H_
#define RUNENGINE_H_

#include<stdio.h>
#include"COORD.h"
#include"PRMTOP.h"
#include <istream>
#include <fstream>

class RunEngine {
public:
	int traj_ans;
	int calc_ans;
	int sel_res;
	int sel_start;
	int sel_end;
	int sel2_start;
	int sel2_end;
	int protein_last_residue;
	bool gzipped;
	bool binnary;
	int sel2_res;
	RunEngine();
	virtual ~RunEngine();
	int Run(char* argv[]);
};

#endif /* RUNENGINE_H_ */
