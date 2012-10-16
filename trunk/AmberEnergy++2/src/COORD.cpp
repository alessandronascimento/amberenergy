/*
 * COORD.cpp
 *
 *  Created on: Dec 21, 2011
 *      Author: asn
 */

#include "COORD.h"

COORD::COORD(PRMTOP* Mol, int as, int ae, int bs, int be, char* filename, bool gzipped) {
	this->astart = as;
	this->aend = ae;
	this->bstart = bs;
	this->bend = be;
	Energy = new ENERGY;
	if (gzipped){
		this->read_gzcrd(Mol, filename);
	}
	else {
		this->read_crd(Mol, filename);
	}
}

COORD::COORD(PRMTOP* Mol, int as, int ae, int bs, int be, char* filename, int mode) {
	this->astart = as;
	this->aend = ae;
	this->bstart = bs;
	this->bend = be;
	Energy = new ENERGY;
	switch (mode){
	case 3:
		this->read_netcdf(Mol, filename);
		break;
	case 4:
		this->read_dcd(Mol, filename);
		break;
	}
}


void COORD::read_crd(PRMTOP* Mol, char* filename) {
	step = 0;
	ifstream mdcrd(filename);
	getline(mdcrd, line);
	printf("# %s\n", line.c_str());
	printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");

	while (!mdcrd.eof()){
		current_crd.clear();
		for (int i=0; i < Mol->N; i++) {
			xyz.clear();
			for (int j=0; j < 3; j++) {
				mdcrd >> pos;
				xyz.push_back(pos);
			}
			current_crd.push_back(xyz);
		}
		if (Mol->ifbox > 0){
			for (int i=0; i<3; i++){
				mdcrd >> pos;	// BOX INFO
			}
		}
		step++;
		printf(" %12d ", step);
		Energy->compute_nb2(Mol, current_crd, this->astart, this->aend, this->bstart, this->bend);
	}
}

void COORD::read_gzcrd(PRMTOP* Mol, char* filename) {
	step = 0;
	igzstream mdcrd(filename);
	getline(mdcrd, line);
	printf("# %s\n", line.c_str());
	printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");

	while (!mdcrd.eof()){
		current_crd.clear();
		for (int i=0; i < Mol->N; i++) {
			xyz.clear();
			for (int j=0; j < 3; j++) {
				mdcrd >> pos;
				xyz.push_back(pos);
			}
			current_crd.push_back(xyz);
		}
		if (Mol->ifbox > 0 ){
			for (int i=0; i< 3; i++){
				mdcrd >> pos; 	// Box Info
			}
		}
		step++;
		printf(" %12d ", step);
		Energy->compute_nb2(Mol, current_crd, this->astart, this->aend, this->bstart, this->bend);
	}
}

COORD::~COORD() {
	current_crd.clear();
	xyz.clear();
}

void COORD::read_netcdf(PRMTOP* Mol, char* filename){
	NcFile nc_mdcrd(filename, NcFile::ReadOnly);

	if (!nc_mdcrd.is_valid()){
		printf("$ Could not open trajectory file %s. Please check.\n", filename);
	}

	NcDim* FrameDim = nc_mdcrd.get_dim("frame");
	int size = FrameDim->size();
	printf("# NetCDF Frame dimension: %d\n", size);

	NcDim* NDim = nc_mdcrd.get_dim("atom");
	const int Ndim = NDim->size();
	if (Ndim != Mol->N){
		printf("# Mismatch among number of atoms in PRMTOP (%d) and NETCDF (%d) files. Please check.\n", Mol->N, Ndim);
		exit(1);
	}
	else {
		printf("# NetCDF number of atoms: %d\n", Ndim);
	}

	NcVar* nc_Coordinates = nc_mdcrd.get_var("coordinates");
	double coords[Ndim][3];

	printf("#%12s %12s %12s %12s\n", "Step", "Elec", "VDW", "Total");

	for (int frame=1; frame <= size; frame++){
		nc_Coordinates->get(&coords[0][0], 1, Ndim, 3);
		printf(" %12d ", frame);
		Energy->compute_nb2(Mol, coords, this->astart, this->aend, this->bstart, this->bend);
		nc_Coordinates->set_cur(frame);
	}
}

void COORD::read_dcd(PRMTOP* Mol, char* filename){
	FILE *dcdfile;
	dcdfile = fopen(filename, "rb");
	int* hdr;
	hdr = new int[1];
	char magic[4];
	int inctrl[20];
	int NTITL;
	int NATREC;
	int NFREAT;
	char* title;
	float *x, *y, *z;
	double *box = new double[6];



	this->parse_binary_int(hdr, dcdfile);

	fread(magic, sizeof(char), 4, dcdfile);
	printf("# DCD HDR: %s\n", magic);

	fread(&inctrl, sizeof(int), 20, dcdfile);

	this->parse_binary_int(hdr, dcdfile);

	this->parse_binary_int(hdr, dcdfile);

	fread(&NTITL, sizeof(int), 1, dcdfile);
	printf("# DCD NTITL: %d\n", NTITL);

	for (int j=0; j < NTITL; j++) {
		title = (char *) malloc(sizeof(char) * 82);
		for (int i=0; i<80; i++) {
			title[i] = fgetc(dcdfile);
		}
		title[80] = '\n';
		title[81] = (char) 0;
		printf("# DCD TITLE: %s\n", title);
	}

	this->parse_binary_int(hdr, dcdfile);

	this->parse_binary_int(hdr, dcdfile);

	fread(&NATREC, sizeof(int), 1, dcdfile);
	printf("# DCD NATREC: %d\n", NATREC);

	this->parse_binary_int(hdr, dcdfile);

	NFREAT = NATREC-inctrl[8];
	printf("# DCD NFREAT: %d\n", NFREAT);

	x = new float[NATREC];
	y = new float[NATREC];
	z = new float[NATREC];

	step = 0;

	while(!feof(dcdfile)){
		step++;

		if (inctrl[10] == 1){
			this->parse_binary_int(hdr, dcdfile);
			if (hdr[0] != -1){
				fread(box, sizeof(double), 6, dcdfile);
//				printf("# DCD Box info: %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", box[0], box[1], box[2], box[3], box[4], box[5]);
				this->parse_binary_int(hdr, dcdfile);
			}
		}

		this->parse_binary_int(hdr, dcdfile);

		fread(x, sizeof(float), NATREC, dcdfile);

		this->parse_binary_int(hdr, dcdfile);
		this->parse_binary_int(hdr, dcdfile);

		fread(y, sizeof(float), NATREC, dcdfile);

		this->parse_binary_int(hdr, dcdfile);
		this->parse_binary_int(hdr, dcdfile);

		fread(z, sizeof(float), NATREC, dcdfile);

		this->parse_binary_int(hdr, dcdfile);

/*
 * Read coordinates. Now, organizing in c++ vector or double array. vector is prefereed.
 */
		xyz.clear();
		current_crd.clear();
		for (int i=0; i< NATREC; i++){
			xyz.clear();
			xyz.push_back(x[i]);
			xyz.push_back(y[i]);
			xyz.push_back(z[i]);
			current_crd.push_back(xyz);
			xyz.clear();
		}
		printf(" %12d ", step);
		Energy->compute_nb2(Mol, current_crd, this->astart, this->aend, this->bstart, this->bend);
	}
	delete x;
	delete y;
	delete z;
	delete box;
}

void COORD::parse_binary_int(int* i, FILE* file){

	fread(i, sizeof(int), 1, file);

#ifdef DEBUG
	printf("HDR: %d\n", i);
#endif

}
