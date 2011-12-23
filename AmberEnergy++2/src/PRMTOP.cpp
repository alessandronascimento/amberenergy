/*
 * PRMTOP.cpp
 *
 *  Created on: 23/03/2010
 *      Author: Alessandro
 */

#include "PRMTOP.h"



PRMTOP::PRMTOP(ifstream &prmtop){
	PRMTOP::get_N(prmtop);
	PRMTOP::get_atomnames(prmtop, PRMTOP::N);
	PRMTOP::get_charges(prmtop, PRMTOP::N);
	PRMTOP::get_masses(prmtop, PRMTOP::N); // unused
	PRMTOP::get_atom_type_index(prmtop, PRMTOP::N);
	PRMTOP::get_nonbond_parm_index(prmtop, PRMTOP::Natomtypes);
	PRMTOP::get_resnames(prmtop);
	PRMTOP::get_res_pointers(prmtop);
	PRMTOP::get_LJA(prmtop, PRMTOP::Natomtypes);
	PRMTOP::get_LJB(prmtop, PRMTOP::Natomtypes);
	PRMTOP::get_H_ACOEF(prmtop, PRMTOP::nph);
	PRMTOP::get_H_BCOEF(prmtop, PRMTOP::nph);
//	PRMTOP::get_amberatoms(prmtop, PRMTOP::N);
//	PRMTOP::read_atomtypes_prm();
//	PRMTOP::get_vdw_parms(PRMTOP::amberatoms, PRMTOP::a_atomtypes_prm, PRMTOP::a_radius, PRMTOP::a_welldepth);
//	PRMTOP::get_coordinates(inpcrd);
}

void PRMTOP::get_N(ifstream &prmtop) {
	if (prmtop.is_open()){

		getline (prmtop, line);

		for (int i=1; i<=5; i++) {
			getline (prmtop, line);
		}

		prmtop >> PRMTOP::N >> PRMTOP::Natomtypes;

		for (int i=1; i<=9; i++){
			prmtop >> line;
		}

		prmtop >> PRMTOP::Nres;

		//
		for (int i=1; i<=7; i++) {
			prmtop >> PRMTOP::nph;
		}

		prmtop >> nph;
		printf("# N: %d Ntypes: %d Nres: %d NPH: %d\n", N, Natomtypes, Nres, nph);

		for (int i=0; i<7; i++){
			prmtop >> this->ifbox;
		}

		prmtop >> ifbox;

		switch (ifbox){
			case 0:
				printf("# There is no periodic box. IFBOX = %d.\n", ifbox);
				break;
			case 1:
				printf("# There is a cubic periodic box. IFBOX = %d.\n", ifbox);
				break;
			case 2:
				printf("# There is a truncated octaedron periodic box. IFBOX = %d.\n", ifbox);
				break;
		}

		getline(prmtop,line);
	}
	else {
		cout << "Couldn't open prmtop file;" << endl;
		exit(1);
	}
}

void PRMTOP::get_atomnames(ifstream &prmtop, int N){
	string name, name2;
	getline (prmtop, line);
	while (line.size() < 6 or line.substr(6,9) != "ATOM_NAME")  {
		getline(prmtop, line);
	}
	getline(prmtop, line);
//	printf("# %s\n", line.c_str());
	int i=0;
	while (i<N){
		prmtop >> name;
		if (name.size() > 4){
			while (name.size()>4) {
				name2=name.substr(0,4);
				PRMTOP::atomnames.push_back(name2);
				i++;
				name=name.substr(4);
			}
			PRMTOP::atomnames.push_back(name);
			i++;
		}
		else {
			PRMTOP::atomnames.push_back(name);
			i++;
		}
	}
}

void PRMTOP::get_charges(ifstream &prmtop, int N) {
	getline (prmtop, line);
	while (line.size() < 6 or line.substr(6,6) != "CHARGE")  {
		getline (prmtop, line);
	}
	getline (prmtop, line);		// FORMAT
//	printf("# %s\n", line.c_str());
	for (int i=1; i<=N; i++) {
		prmtop >> charge;
		PRMTOP::charges.push_back(charge);
	}

	getline (prmtop, line);
}


void PRMTOP::get_masses(ifstream &prmtop, int N) {
	getline (prmtop, line);
	while (line.size() < 6 or line.substr(6,4) != "MASS")  {
		getline (prmtop, line);
	}
	getline (prmtop, line);		// FORMAT
//	printf("# %s\n", line.c_str());
	for (int i=1; i<=N; i++) {
		prmtop >> mass;
		PRMTOP::masses.push_back(mass);
	}
	getline (prmtop, line);
}

void PRMTOP::get_resnames(ifstream &prmtop){
	getline (prmtop, line);
	string resname;
	while (line.size() < 6 or line.substr(6,13) != "RESIDUE_LABEL")  {
		getline (prmtop, line);
	}
	getline (prmtop, line);		// FORMAT
//	printf("# %s\n", line.c_str());
	for (int i=1; i<= this->Nres; i++){
		prmtop >> resname;
		this->resnames.push_back(resname);
	}
	getline (prmtop, line);
}

void PRMTOP::get_res_pointers(ifstream &prmtop){
	getline(prmtop, line);
	int respointer;
	while (line.size() < 6 or line.substr(6,15)!= "RESIDUE_POINTER")  {
		getline (prmtop, line);
	}
	getline (prmtop, line);		// FORMAT
//	printf("# %s\n", line.c_str());
	for (int i=0; i < this->Nres; i++){
		prmtop >> respointer;
		PRMTOP::residue_pointer.push_back(respointer);
	}
	getline(prmtop, line);
}

void PRMTOP::get_atom_type_index(ifstream &prmtop, int N){

	int atomtype;

	getline(prmtop, line);

	while (line.substr(6,15) != "ATOM_TYPE_INDEX")  {
		getline(prmtop, line);
	}

	getline (prmtop, line);		// FORMAT
//	printf("# %s\n", line.c_str());
		//  cout << line << endl;

	for (int i=1; i<=N; i++) {
		prmtop >> atomtype;
		this->atom_types_index.push_back(atomtype);
	}
	getline (prmtop, line);
}

void PRMTOP::get_nonbond_parm_index(ifstream &prmtop, int Natomtypes){
	getline (prmtop, line);
	int nb_parm_index;
	while (line.substr(6,20) != "NONBONDED_PARM_INDEX") {
		getline(prmtop, line);
	}
	getline(prmtop, line);		// FORMAT
//	printf("# %s\n", line.c_str());

	for (int i=1; i<= (Natomtypes*Natomtypes); i++) {
			prmtop >> nb_parm_index;
			this->ico.push_back(nb_parm_index);
	}
		getline (prmtop, line);
}

void PRMTOP::get_amberatoms(ifstream &prmtop, int N) {
	getline (prmtop, line);
	while (! prmtop.eof()) {
		getline (prmtop, line);
		if (line.size() >= 15) {
			if (line.substr(6,15) == "AMBER_ATOM_TYPE")  {
				getline (prmtop, line);		// FORMAT
//				printf("# %s\n", line.c_str());
				for (int i=1; i<=N; i++) {
					prmtop >> atm;
					PRMTOP::amberatoms.push_back(atm); }
			}
		}
	}
}

void PRMTOP::get_LJA(ifstream &prmtop, int Natomtypes){
	double LJ_a;
	getline (prmtop, line);
	while (line.substr(0,5) != "%FLAG" || line.substr(6,19) != "LENNARD_JONES_ACOEF") {
		getline(prmtop, line);
	}

	getline(prmtop, line);		//FORMAT
//	printf("# %s\n", line.c_str());

	for (int i=1; i<= (Natomtypes*(Natomtypes+1)/2); i++) {
		prmtop >> LJ_a;
		this->CN1.push_back(LJ_a);
	}
		getline(prmtop, line);
}

void PRMTOP::get_LJB(ifstream &prmtop, int Natomtypes){
	double LJ_b;
	getline(prmtop, line);
	while (line.substr(0,5) != "%FLAG" || line.substr(6,19) != "LENNARD_JONES_BCOEF") {
		getline(prmtop, line);
	}

	getline(prmtop, line);		//FORMAT
//	printf("# %s\n", line.c_str());

	for (int i=1; i<= (Natomtypes*(Natomtypes+1)/2); i++) {
		prmtop >> LJ_b;
		this->CN2.push_back(LJ_b);
	}
		getline(prmtop, line);
}

void PRMTOP::get_H_ACOEF(ifstream &prmtop, int nph){
	double asol;
	getline(prmtop, line);
	while (line.substr(0,5) != "%FLAG" || line.substr(6,11) != "HBOND_ACOEF") {
		getline(prmtop, line);
	}

	getline(prmtop, line);		//FORMAT
//	printf("# %s\n", line.c_str());

	for (int i=1; i<= nph; i++) {
		prmtop >> asol;
		this->ASOL.push_back(asol);
	}

	getline(prmtop, line);
}

void PRMTOP::get_H_BCOEF(ifstream &prmtop, int nph){
	double bsol;
	while (line.substr(0,5) != "%FLAG" || line.substr(6,11) != "HBOND_BCOEF") {
		getline(prmtop, line);
	}

	getline(prmtop, line);		//FORMAT
//	printf("# %s\n", line.c_str());

	for (int i=1; i<= nph; i++) {
		prmtop >> bsol;
		BSOL.push_back(bsol);
	}

	getline(prmtop, line);
}
