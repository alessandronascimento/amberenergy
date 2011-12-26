/*
 * PRMTOP.h
 *
 *  Created on: 23/03/2010
 *      Author: Alessandro
 */

#ifndef PRMTOP_H_
#define PRMTOP_H_

#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<stdio.h>
#include<stdlib.h>

using namespace std;
/*!
 * The PRMTOP class has a number of methods dedicated to parse a AMBER prmtop
 * file and get all the necessary information for energy evaluation in a receptor-ligand
 * system. There is also a method to parse the atomic coordinates of a inpcrd file. The
 * parser works well for ligand and receptor prmtop/inpcrd files.
 */
class PRMTOP {

public:

	// variables

	//! Number of atoms
	int N;
	//! Number of residues
	int Nres;
	//! Number of AMBER atomtypes
	int Natomtypes;
	//! number of excluded atoms
	int nnb;
	//! number of distinct 10-12 hydrogen bond pair types
	int nph;
	//!
	int ifbox;
	//! Counter variable
	int count;
	//! String to parse atom names
	string atm,
	//! String to parse comment lines
	line;
	//!
	vector<string> atomnames;
	//! C++ vector to store atomic charges
	vector<double>charges;
	//! C++ vector to store atomic masses
	vector<double>masses;
	//! C++ vector to store the residue labels
	vector<string> resnames;
	//!
	vector<int> atom_types_index;
	//!
	vector<int>number_excluded_atoms;
	//!
	vector<int> excluded_atoms_list;
	//!
	vector<int> ico;
	//! C++ vector to store the AMBERatoms
	vector<string>amberatoms;
	//! C++ vector with pointer to the number of the atoms in the residue
	vector<int> residue_pointer;
	//! Atomic mass
	double mass,
	//! Atomic charge
	charge;
	//!
	vector<double> CN1;
	//!
	vector<double> CN2;
	//!
	vector<double> ASOL;
	//!
	vector<double> BSOL;


	//methods used to parse prmtop/inpcrd file
/*!
 * Class constructor. The some internal methods are called to parse the
 * number of atoms, the atomic charges, the atomic masses, the amberatom
 * names, and the coordinates. Also, the parameters of the force field are
 * taken from the file "vdw.param" (epsilons, radius)
 */
	PRMTOP(ifstream &prmtop);

	/*!
	 * This method parses the prmtop file to get the number
	 * of atoms in the system
	 */
	void get_N(ifstream &prmtop);

	//!
	void get_atomnames(ifstream &prmtop, int N);

	/*!
	 * This method parses the prmtop file to get the atomic
	 * charges (electron charge multiplied by 18.2223).
	 */
	void get_charges(ifstream &prmtop, int N);

	/*!
	 * This method parses the PRMTOP file to get the residue
	 * label.
	 */
	void get_resnames(ifstream &prmtop);

	/*!
	 * This method parses the PRMTOP file to get the pointers to
	 * residues. That is, the number of the initial atom of each residue.
	 */
	void get_res_pointers(ifstream &prmtop);

	/*!
	 * This method parses the prmtop file to get the atomic
	 * masses of the atoms in the system and stores them in a
	 * C++ vector
	 */
	void get_masses(ifstream &prmtop, int N);

	/*!
	 * This method parses the atomic types of the system in
	 * Amber FF type.
	 */
	void get_amberatoms(ifstream &prmtop, int N);


	void get_atom_type_index(ifstream &prmtop, int N);

	void get_nonbond_parm_index(ifstream &prmtop, int Natomtypes);

	void get_LJA(ifstream &prmtop, int Natomtypes);

	void get_LJB(ifstream &prmtop, int Natomtypes);

	void get_H_ACOEF(ifstream &prmtop, int nph);

	void get_H_BCOEF(ifstream &prmtop, int nph);

	void get_number_excluded_atoms(ifstream &prmtop, int N);
	void get_excluded_atoms_list(ifstream &prmtop, int nnb);
};

#endif /* PRMTOP_H_ */
