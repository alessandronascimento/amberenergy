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
#include<vector>
#include<math.h>

using namespace std;

double x, y, z, r, r2, elec, vdw;
double diel=1.00;
int iaci, ic;


double compute_r(vector<double> i, vector<double> j) {
	x= j[0]-i[0];
	y= j[1]-i[1];
	z= j[2]-i[2];
	r= sqrt((x*x)+(y*y)+(z*z)) ;
	return (r) ; }

double compute_nb (vector<vector<double> > current, int astart, int aend, vector<int> iac, vector<int> ico) {
	elec=0.00;
	vdw=0.00;
	for (int i= astart; i<= aend; i++) {
		for (int j=0; j<current.size(); j++) {
			if (j < astart || j > aend) {
				r = compute_r(current[i], current[j]);
				r2 = r*r;
				elec = elec+((charges[i]*charges[j])/(r*diel));
				iaci = Natomtypes*(iac[i]-1);
				ic = ico[(iaci+iac[j])-1];
				if (ic > 0) { 		// Use 12 - 6 LJ Potencial
					vdw = vdw + (((LJA[ic-1])/(r2*r2*r2*r2*r2*r2))-((LJB[ic-1])/(r2*r2*r2))); }
				else {			// Use 10 - 12 LJ Potencial
					vdw = vdw + (((LJA[ic-1])/(r2*r2*r2*r2*r2*r2))-((LJB[ic-1])/(r2*r2*r2*r2*r2))); } } } }
	printf("%12.4f %12.4f %12.4f\n", elec, vdw, (elec+vdw));
	return (elec+vdw);}

double compute_nb2 (vector<vector<double> > current, int astart, int aend, vector<int> iac, vector<int> ico, int bstart, int bend) {
	elec=0.00;
	vdw=0.00;
	for (int i= astart; i<= aend; i++) {
		for (int j=bstart; j<=bend; j++) {
			r = compute_r(current[i], current[j]);
			r2 = r*r;
			elec = elec+((charges[i]*charges[j])/(r*diel));
			iaci = Natomtypes*(iac[i]-1);
			ic = 	ico[(iaci+iac[j])-1];
			if (ic > 0) { 		// Use 12 - 6 LJ Potencial
				vdw = vdw + (((LJA[ic-1])/(r2*r2*r2*r2*r2*r2))-((LJB[ic-1])/(r2*r2*r2))); }

			else {			// Use 10 - 12 LJ Potencial
				vdw = vdw + (((ASOL[-(ic-1)])/(r2*r2*r2*r2*r2*r2)) - ((BSOL[-(ic-1)])/(r2*r2*r2*r2*r2)));} } }
	printf("%12.4f %12.4f %12.4f\n", elec, vdw, (elec+vdw));
	return (elec+vdw);}
