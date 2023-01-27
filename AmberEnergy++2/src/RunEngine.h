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
#include<vector>

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
    int receptor_chuncks;
    int ligand_chuncks;
//    vector<vector<int> > receptor_residues_chuncks;
//    vector<vector<int> > receptor_atom_chuncks;
//    vector<vector<int> > ligand_residues_chuncks;
//    vector<vector<int> > ligand_atom_chuncks;
    vector<int> receptor_atoms;
    vector<int> ligand_atoms;
	RunEngine();
	virtual ~RunEngine();
	int Run(char* argv[]);
};

#endif /* RUNENGINE_H_ */
