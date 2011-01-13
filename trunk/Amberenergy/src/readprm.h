/*
 * readprm.cpp
 *
 *  Created on: 16/09/2010
 *      Author: Nascimento
 */

#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;


int N, Natomtypes, atomtype, nb_parm_index, Nres, residt, nph;
double charge, LJ_a, LJ_b, asol, bsol;
string line, flag, resn;
vector<double> charges, LJA, LJB, ASOL, BSOL;
vector<int> atom_types, ico, res_pointer;
vector<string> resnames;


void read_parm(ifstream &prmtop) {
	getline (prmtop, line);		// 1st line
	cout << "# " << line << endl;
	for (int i=1; i<=5; i++) {
		getline (prmtop, line); }
	prmtop >> N >> Natomtypes;
	for (int i=1; i<=9; i++) {
		prmtop >> Nres; }
	prmtop >> Nres;
	for (int i=1; i<=7; i++) {
		prmtop >> nph; }
	prmtop >> nph;
	printf("# N: %d Ntypes: %d Nres: %d NPH: %d\n", N, Natomtypes, Nres, nph);
	getline(prmtop,line);
	getline(prmtop,line);
	/*

Getting the charges using amber default format

	 */
	while (line.substr(6,6) != "CHARGE")  {
		getline (prmtop, line); }
	getline (prmtop, line);		// FORMAT
	//  cout << line << endl;
	for (int i=1; i<=N; i++) {
		prmtop >> charge;
		charges.push_back(charge); }
	getline (prmtop, line);
	getline (prmtop, line);

	/*

Getting the atomtypes indexes using amber default format

	 */

	while (line.substr(6,15) != "ATOM_TYPE_INDEX")  {
		getline(prmtop, line); }
	getline (prmtop, line);		// FORMAT
	//  cout << line << endl;
	for (int i=1; i<=N; i++) {
		prmtop >> atomtype;
		atom_types.push_back(atomtype); }
	getline (prmtop, line);
	getline (prmtop, line);

	/*

Getting the ico indexes using amber default format

	 */

	while (line.substr(6,20) != "NONBONDED_PARM_INDEX") {
		getline(prmtop, line); }
	getline(prmtop, line);		// FORMAT
	//  cout << line << endl;
	for (int i=1; i<= (Natomtypes*Natomtypes); i++) {
		prmtop >> nb_parm_index;
		ico.push_back(nb_parm_index); }
	getline (prmtop, line);
	getline (prmtop, line);

	/*

Getting the residue names

	 */
	while (line.substr(6,13) != "RESIDUE_LABEL") {
		getline(prmtop, line); }
	getline(prmtop, line);		// FORMAT
	//	cout << line << endl;
	for (int i=1; i<= Nres; i++) {
		prmtop >> resn;
		resnames.push_back(resn); }
	getline(prmtop, line);
	getline(prmtop, line);


	/*

Getting the residue pointer indexes using amber default format

	 */

	while (line.substr(6,15) != "RESIDUE_POINTER") {
		getline(prmtop, line); }
	getline(prmtop, line);		// FORMAT
	//  cout << line << endl;
	for (int i=1; i<= Nres; i++) {
		prmtop >> residt;
		res_pointer.push_back(residt); }
	getline(prmtop, line);
	getline(prmtop, line);
	/*

Getting the LJ indexes using amber default format

	 */
	while (line.substr(0,5) != "%FLAG" || line.substr(6,19) != "LENNARD_JONES_ACOEF") {
		getline(prmtop, line); }
	getline(prmtop, line);		//FORMAT
	//  cout << line << endl;
	for (int i=1; i<= (Natomtypes*(Natomtypes+1)/2); i++) {
		prmtop >> LJ_a;
		LJA.push_back(LJ_a); }
	getline(prmtop, line);
	getline(prmtop, line);

	while (line.substr(0,5) != "%FLAG" || line.substr(6,19) != "LENNARD_JONES_BCOEF") {
		getline(prmtop, line); }
	getline(prmtop, line);		//FORMAT
	//  cout << line << endl;
	for (int i=1; i<= (Natomtypes*(Natomtypes+1)/2); i++) {
		prmtop >> LJ_b;
		LJB.push_back(LJ_b); }
	getline(prmtop, line);
	getline(prmtop, line);


	while (line.substr(0,5) != "%FLAG" || line.substr(6,11) != "HBOND_ACOEF") {
		getline(prmtop, line); }
	getline(prmtop, line);		//FORMAT
	for (int i=1; i<= nph; i++) {
		prmtop >> asol;
		ASOL.push_back(asol); }
	getline(prmtop, line);
	getline(prmtop, line);

	while (line.substr(0,5) != "%FLAG" || line.substr(6,11) != "HBOND_BCOEF") {
		getline(prmtop, line); }
	getline(prmtop, line);		//FORMAT
	for (int i=1; i<= nph; i++) {
		prmtop >> bsol;
		BSOL.push_back(bsol); }
	getline(prmtop, line);
	getline(prmtop, line);

}
