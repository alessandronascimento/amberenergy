/*
 * readcrd.cpp
 *
 *  Created on: 16/09/2010
 *      Author: Nascimento
 */

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
#include<fstream>
#include<string>
#include<vector>
#include<stdio.h>
#include "gzstream.h"
#include "process_nb.h"

using namespace std;

int step=0;
double pos;
vector<vector<double> >current_crd;
vector<double> xyz;

void read_crd(int sel_start, int sel_end, char* argv[]) {
  ifstream mdcrd(argv[2]);
  getline(mdcrd, line);
  cout << "# " << line << endl;
  printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");
  while (!mdcrd.eof()){
    current_crd.clear();
    for (int i=1; i<=N; i++) {
      xyz.clear();
      for (int j=1; j<=3; j++) {
	mdcrd >> pos;
	xyz.push_back(pos); }
      current_crd.push_back(xyz); }
    step++;
    printf(" %12d ", step);
    double nonbonded = compute_nb(current_crd, sel_start-1, sel_end-1, atom_types, ico);  } }

void read_gzcrd(int sel_start, int sel_end, char* argv[]) {
  igzstream mdcrd(argv[2]);
  getline (mdcrd, line);
  cout << "# " << line << endl;
  printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");
  while (!mdcrd.eof()){
    current_crd.clear();
    for (int i=1; i<=N; i++) {
      xyz.clear();
      for (int j=1; j<=3; j++) {
	mdcrd >> pos;
	xyz.push_back(pos); }
	current_crd.push_back(xyz); }
    step++;
    printf(" %12d ", step);
    double elect = compute_nb(current_crd, sel_start-1, sel_end-1, atom_types, ico); } }

void read_crd_box(int sel_start, int sel_end, char* argv[]) {
  ifstream mdcrd(argv[2]);
  getline(mdcrd, line);
  cout << "# " << line << endl;
  printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");
  while (!mdcrd.eof()){
    current_crd.clear();
    for (int i=1; i<=N; i++) {
      xyz.clear();
      for (int j=1; j<=3; j++) {
	mdcrd >> pos;
	xyz.push_back(pos); }
      current_crd.push_back(xyz); }
    for (int i=1; i<=3; i++) {		// BOX INFO
      mdcrd >> pos; }
    step++;
    printf(" %12d ", step);
    double elect = compute_nb(current_crd, sel_start-1, sel_end-1, atom_types, ico); } }

void read_gzcrd_box(int sel_start, int sel_end, char* argv[]) {
  igzstream mdcrd(argv[2]);
  getline (mdcrd, line);
  cout << "# " << line << endl;
  printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");
  while (!mdcrd.eof()){
    current_crd.clear();
    for (int i=1; i<=N; i++) {
      xyz.clear();
      for (int j=1; j<=3; j++) {
	mdcrd >> pos;
	xyz.push_back(pos); }
      current_crd.push_back(xyz); }
    for (int i=1; i<=3; i++) {		// BOX INFO
      mdcrd >> pos; }
    step++;
    printf(" %12d ", step);
    double elect = compute_nb(current_crd, sel_start-1, sel_end-1, atom_types, ico); } }


void read_crd2(int sel_start, int sel_end, int sel2_start, int sel2_end, char* argv[]) {
  ifstream mdcrd(argv[2]);
  getline(mdcrd, line);
  cout << "# " << line << endl;
  printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");
  while (!mdcrd.eof()){
  current_crd.clear();
  for (int i=1; i<=N; i++) {
      xyz.clear();
      for (int j=1; j<=3; j++) {
	mdcrd >> pos;
	xyz.push_back(pos); }
    current_crd.push_back(xyz); }
    step++;
    printf(" %12d ", step);
    double elect = compute_nb2(current_crd, sel_start-1, sel_end-1, atom_types, ico, sel2_start-1, sel2_end-1); } }

void read_gzcrd2(int sel_start, int sel_end, int sel2_start, int sel2_end, char* argv[]) {
  igzstream mdcrd(argv[2]);
  getline (mdcrd, line);
  cout << "# " << line << endl;
  printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");
  while (!mdcrd.eof()){
    current_crd.clear();
    for (int i=1; i<=N; i++) {
      xyz.clear();
      for (int j=1; j<=3; j++) {
	mdcrd >> pos;
	xyz.push_back(pos); }
      current_crd.push_back(xyz); }
    step++;
    printf(" %12d ", step);
    double elect = compute_nb2(current_crd, sel_start-1, sel_end-1, atom_types, ico, sel2_start-1, sel2_end-1); } }

void read_crd_box2(int sel_start, int sel_end, int sel2_start, int sel2_end, char* argv[]) {
  ifstream mdcrd(argv[2]);
  getline(mdcrd, line);
  cout << "# " << line << endl;
  printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");
  while (!mdcrd.eof()){
  current_crd.clear();
  for (int i=1; i<=N; i++) {
      xyz.clear();
      for (int j=1; j<=3; j++) {
	mdcrd >> pos;
	xyz.push_back(pos); }
    current_crd.push_back(xyz); }
    for (int k=1; k<=3; k++) {		// BOX INFO
      mdcrd >> pos; }
    step++;
    printf(" %12d ", step);
    double elect = compute_nb2(current_crd, sel_start-1, sel_end-1, atom_types, ico, sel2_start-1, sel2_end-1); } }

void read_gzcrd_box2(int sel_start, int sel_end, int sel2_start, int sel2_end, char* argv[]) {

	igzstream mdcrd(argv[2]);
	float temp;

	getline (mdcrd, line);
	cout << "# " << line << endl;

	printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");

	while (!mdcrd.eof()){

		current_crd.clear();
		for (int i=0; i<N; i++) { // all atoms in a frame
			for (int j=0; j<3; j++) { // x, y and z coordinates
				mdcrd >> pos;
				xyz.push_back(pos);
			}
			current_crd.push_back(xyz);
			xyz.clear();
		}

		if (ifbox > 0){
//			printf ("Box info: ");
			for (int k=1; k<=3; k++) {		// BOX INFO
				mdcrd >> temp;
//				printf("%.4f ", pos);
			}
//			printf("\n");
		}

		step++;

		printf(" %12d ", step);
		double elect = compute_nb2(current_crd, sel_start-1, sel_end-1, atom_types, ico, sel2_start-1, sel2_end-1);
	}
}
