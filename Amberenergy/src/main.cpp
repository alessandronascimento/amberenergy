/*
 * main.cpp
 *
 *  Created on: 16/09/2010
 *      Author: Nascimento
 */

/*

/*

#------------------------------------------------------------------------#
#                                                                        #
# AMBERENERGY - A PROGRAM TO EVALUATE ENERGIES IN A AMBER MM TRAJECTORY  #
#             Written by Alessandro S. Nascimento -   Aug/2008           #
#                      al.s.nascimento at gmail.com                      #
#                                                                        #
#------------------------------------------------------------------------#

*/

#include<iostream>
#include<stdlib.h>
#include<fstream>
#include<string>
#include<vector>
#include<ctime>
#include "readprm.h"
#include "readcrd.h"
#include <gzstream.h>

using namespace std;

int sel_res, sel_start, sel_end, sel2_res, sel2_start, sel2_end, traj_ans, calc_ans;
//string ans_yn;
clock_t start, finish;

int main (int argc, char* argv[]) {

	if (argc < 2) {
		cout << "USAGE: amberenergy.x prmtop mdcrd" << endl;
		exit(1); }

	cout << "# Reading prmtop information..." << endl;
	ifstream prmtop (argv[1]);
	read_parm(prmtop);
	prmtop.close();
	cout <<"# Prmtop file " << argv[1] << " read!" << endl;

	/*   TRAJECTORY FILE FORMAT    */

	cout << endl;
	cout << "# Define the format of the trajectory file: " << endl;
	cout << endl;
	cout << "# \t 1) AMBER MDCRD FILE WITH NO BOX" << endl;
	cout << "# \t 2) GZIPPED AMBER MDCRD FILE WITH NO BOX" << endl;
	cout << "# \t 3) AMBER MDCRD FILE WITH BOX INFORMATION" << endl;
	cout << "# \t 4) GZIPPED AMBER MDCRD FILE WITH BOX INFORMATION" << endl;
	cout << "# \t 5) NETCDF FORMAT FILE (not implemented yet)" << endl;
	cout << endl;
	cout << "# Enter your choice number: " ;
	cin >> traj_ans;


	/* SELECTIONS */

	cout << endl;
	cout <<"# Choose your option for energy calculations: ";
	cout << endl;
	cout << "# \t 1) Ligand - Protein interaction calculations" << endl;
	cout << "# \t 2) Ligand - Environment interaction calculations" << endl;
	cout << "# \t 3) Ligand - Solvent interaction calculations" << endl;
	cout << "# \t 4) Selection - Selection interaction calculations " << endl;
	cout << "# \t 5) Protein Internal Energy Calculations" << endl;
	cout << "# \t 6) Environment Internal Energy Calculations" << endl;
	cout << "# \t 7) Multiple Selection - Multiple Selection Calculations" << endl;
	cout << endl;
	cout << "# Enter your choice number: " ;
	cin >> calc_ans;

	int protein_last_residue = 0;
	switch(calc_ans) {

	case 2:
		cout << "# Define your ligand residue number: " ;
		cin >> sel_res;
		sel_start = res_pointer[sel_res-1];
		sel_end = res_pointer[sel_res]-1;
		cout << "# You selected atoms " << sel_start << " to " << sel_end << endl;
		start = clock();
		if (traj_ans == 1) {
			read_crd(sel_start, sel_end, argv); }
		else if (traj_ans == 3) {
			read_crd_box(sel_start, sel_end, argv); }
		else if (traj_ans == 2) {
			read_gzcrd(sel_start, sel_end, argv); }
		else if (traj_ans == 4) {
			read_gzcrd_box(sel_start, sel_end, argv); }
		break;

	case 1:
		cout << "# Define your ligand residue number: " ;
		cin >> sel_res;
		sel_start = res_pointer[sel_res-1];
		sel_end = res_pointer[sel_res]-1;
		for(int i=0; i<resnames.size(); i++) {
			if (resnames[i] != "WAT" and resnames[i] != "Na+" and resnames[i] != "Cl-" and resnames[i] != "Mn2" and resnames[i] != "Ca2") {
				protein_last_residue++; } }
		cout << "# Protein last residue number: " << protein_last_residue << endl;
		sel2_start = sel_end+1;
		sel2_end = res_pointer[protein_last_residue]-1;
		if (sel_end > sel2_end) {				// protein comes before ligand
			sel2_start = 1; }
		cout << "# You selected atom " << sel2_start << "(" << atom_types[sel2_start-1]<< ") to " << sel2_end << "(" <<  atom_types[sel2_end-1] << ")" << endl;
		cout << "# And atom " << sel_start << "(" <<  atom_types[sel_start-1] << ") to " << sel_end << "(" <<  atom_types[sel_end-1] << ")" << endl;
		if (traj_ans == 1) {
			read_crd2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 3) {
			read_crd_box2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 2) {
			read_gzcrd2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 4) {
			read_gzcrd_box2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		break;

	case 3:
		cout << "# Define your ligand residue number: " ;
		cin >> sel_res;
		sel_start = res_pointer[sel_res-1];
		sel_end = res_pointer[sel_res]-1;
		cout << "# You selected atoms " << sel_start << " to " << sel_end << endl;
		for(int i=0; i<resnames.size(); i++) {
			if (resnames[i] != "WAT" and resnames[i] != "Na+" and resnames[i] !="Cl-") {
				protein_last_residue++; } }
		cout << "# Protein last residue number : " << protein_last_residue << endl;
		sel2_start = res_pointer[protein_last_residue]; 	// first solvent atom
		sel2_end = N;
		cout << "# You selected atom " << sel2_start << " to " << sel2_end << endl;
		if (traj_ans == 1) {
			read_crd2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 3) {
			read_crd_box2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 2) {
			read_gzcrd2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 4) {
			read_gzcrd_box2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		break;

	case 4:
		cout << "# Define your ligand residue number: " ;
		cin >> sel_res;
		sel_start = res_pointer[sel_res-1];
		sel_end = res_pointer[sel_res]-1;
		cout << "# You selected atom " << sel_start << " to " << sel_end << endl;
		cout << "# Define your second residue: " ;
		cin >> sel2_res;
		sel2_start = res_pointer[sel2_res-1];
		sel2_end = res_pointer[sel2_res]-1;
		cout << "# You selected atom " << sel2_start << " to " << sel2_end << endl;
		if (traj_ans == 1) {
			read_crd2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 3) {
			read_crd_box2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 2) {
			read_gzcrd2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 4) {
			read_gzcrd_box2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		break;

	case 5:
		for(int i=0; i<resnames.size(); i++) {
			if (resnames[i] != "WAT" and resnames[i] != "Na+" and resnames[i] !="Cl-") {
				protein_last_residue++; } }
		cout << "# Protein last residue number : " << protein_last_residue << endl;
		sel2_start = res_pointer[0];			// first solvent atom
		sel2_end = res_pointer[protein_last_residue];
		cout << "# You selected atom " << sel2_start << " to " << sel2_end << endl;
		if (traj_ans == 1) {
			read_crd(sel2_start, sel2_end, argv); }
		else if (traj_ans == 3) {
			read_crd_box(sel2_start, sel2_end, argv); }
		else if (traj_ans == 2) {
			read_gzcrd(sel2_start, sel2_end, argv); }
		else if (traj_ans == 4) {
			read_gzcrd_box(sel2_start, sel2_end,argv); }
		break;

	case 6:
		sel2_start = res_pointer[0];			// first solvent atom
		sel2_end = N;
		cout << "# You selected atom " << sel2_start << " to " << sel2_end << endl;
		if (traj_ans == 1) {
			read_crd(sel2_start, sel2_end, argv); }
		else if (traj_ans == 3) {
			read_crd_box(sel2_start, sel2_end, argv); }
		else if (traj_ans == 2) {
			read_gzcrd(sel2_start, sel2_end, argv); }
		else if (traj_ans == 4) {
			read_gzcrd_box(sel2_start, sel2_end,argv); }
		break;

	case 7:
		cout << "# Define your first residue range (e.g: 1 100): " ;
		cin >> sel_res >> sel2_res ;
		sel_start = res_pointer[sel_res-1];
		sel_end = res_pointer[sel2_res]-1;
		cout << "# You selected atom " << sel_start << " to " << sel_end << endl;
		cout << "# Define your second residue range (e.g: 101 200): " ;
		cin >> sel_res >> sel2_res ;
		sel2_start = res_pointer[sel_res-1];
		sel2_end = res_pointer[sel2_res]-1;
		cout << "# You selected atom " << sel2_start << " to " << sel2_end << endl;
		if (traj_ans == 1) {
			read_crd2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 3) {
			read_crd_box2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 2) {
			read_gzcrd2(sel_start, sel_end, sel2_start, sel2_end, argv); }
		else if (traj_ans == 4) {
			read_gzcrd_box2(sel_start, sel_end, sel2_start, sel2_end, argv); }

		break; }
	finish = clock();
	cout << "# Time elapsed for trajectory analysis: " << (double(finish)-double(start))/CLOCKS_PER_SEC << " seconds." << endl;
	cout << "# Average time/frame: " << ((double(finish)-double(start))/CLOCKS_PER_SEC)/step << " seconds." << endl;
}
