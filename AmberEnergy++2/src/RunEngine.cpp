/*
 * RunEngine.cpp
 *
 *  Created on: Dec 21, 2011
 *      Author: asn
 */

#include "RunEngine.h"

RunEngine::RunEngine() {
    this->gzipped = false;
    this->binnary = false;
}

RunEngine::~RunEngine() {
    // TODO Auto-generated destructor stub
}

int RunEngine::Run(char* argv[]){
    printf("# Reading prmtop information...\n");

    ifstream prmtop (argv[1]);
    if (!prmtop.is_open()){
        printf("# Prmtop file %s could not be opened. Please check.\n", argv[1]);
        exit(1);
    }

    PRMTOP* Mol = new PRMTOP(prmtop);
    printf("# Prmtop file %s read successfully!\n", argv[1]);

    prmtop.close();



    /*   TRAJECTORY FILE FORMAT    */

    printf("#\n");
    printf("# Define the format of the trajectory file: \n");
    printf("#\n");
    printf("# \t 1) AMBER MDCRD FILE\n");
    printf("# \t 2) GZIPPED AMBER MDCRD FILE\n");
    printf("# \t 3) NETCDF FORMAT FILE\n");
    printf("# \t 4) CHARMM DCD FORMAT FILE\n");
    printf("#\n");
    printf("# Enter your choice number: ");
    cin >> traj_ans;

    if (traj_ans <1 or traj_ans > 4){
        printf ("Selection %d invalid. Please try again.\n", traj_ans);
        exit(1);
    }
    if (traj_ans == 2){
        gzipped=true;
    }
    else if (traj_ans == 3){
        this->binnary = true;
    }
    else if (traj_ans == 4){
        this->binnary = true;
    }



    /* SELECTIONS */

    printf("#\n");
    printf("# Choose your option for energy calculations:\n");
    printf("#\n");
    printf("# \t 1) Ligand - Protein interaction calculations\n");
    printf("# \t 2) Ligand - Environment interaction calculations\n");
    printf("# \t 3) Ligand - Solvent interaction calculations\n");
    printf("# \t 4) Selection - Selection interaction calculations\n");
    printf("# \t 5) Protein Internal Energy Calculations\n");
    printf("# \t 6) Multiple Selection - Multiple Selection Calculations\n");
    printf("# \t 7) Multiple Selection - Multiple Selection w/ Chuncks of Residues for Calculations\n");
    printf("#\n");
    printf("# Enter your choice number: ");
    cin >> calc_ans;

    protein_last_residue = 0;

    COORD* Coord;
    vector<int> vitmp;


    switch(calc_ans) {

    case 1:
        cout << "# Define your ligand residue number: " ;
        cin >> sel_res;

        sel_start = Mol->residue_pointer[sel_res-1];
        sel_end = Mol->residue_pointer[sel_res]-1;

        for(unsigned i=0; i<Mol->resnames.size(); i++) {
            if (Mol->resnames[i] != "WAT" and Mol->resnames[i] != "Na+" and Mol->resnames[i] != "Cl-" and Mol->resnames[i] != "Mn2" and Mol->resnames[i] != "Ca2") {
                protein_last_residue++;
            }
        }
        printf("# Protein last residue number: %d\n", protein_last_residue);

        sel2_start = sel_end+1;

        sel2_end = Mol->residue_pointer[protein_last_residue]-1;

        if (sel_end > sel2_end) {				// protein comes before ligand
            sel2_start = 1;
        }

        printf("# You selected atom %d (%s) to %d (%s)\n", sel2_start, Mol->atomnames[sel2_start-1].c_str(), sel2_end, Mol->atomnames[sel2_end-1].c_str());
        printf("# And atom %d (%s) to %d (%s).\n", sel_start, Mol->atomnames[sel_start-1].c_str(), sel_end, Mol->atomnames[sel_end-1].c_str());
        if (this->binnary){
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], traj_ans);
        }
        else {
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], gzipped);
        }

        break;


    case 2:
        printf("# Define your ligand residue number: ");
        cin >> sel_res;
        sel_start = Mol->residue_pointer[sel_res-1];
        sel_end = Mol->residue_pointer[sel_res]-1;
        printf("# You selected atoms %d (%s) to %d (%s)\n", sel_start, Mol->atomnames[sel_start-1].c_str(), sel_end, Mol->atomnames[sel_end-1].c_str());
        sel2_start = sel_end +1;
        sel2_end = Mol->N;
        printf("# And atom %d (%s) to %d (%s).\n", sel_start, Mol->atomnames[sel_start-1].c_str(), sel_end, Mol->atomnames[sel_end-1].c_str());
        if (this->binnary){
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], traj_ans);
        }
        else {
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], gzipped);
        }

        break;



    case 3:
        printf("# Define your ligand residue number: ");
        cin >> sel_res;
        sel_start = Mol->residue_pointer[sel_res-1];
        sel_end = Mol->residue_pointer[sel_res]-1;
        printf("# You selected atoms %d (%s) to %d (%s)\n", sel_start, Mol->atomnames[sel_start-1].c_str(), sel_end, Mol->atomnames[sel_end-1].c_str());

        for(unsigned i=0; i<Mol->resnames.size(); i++) {

            if (Mol->resnames[i] != "WAT" and Mol->resnames[i] != "Na+" and Mol->resnames[i] !="Cl-" and Mol->resnames[i] != "Mn2" and Mol->resnames[i] != "Ca2") {
                protein_last_residue++;
            }
        }

        printf("# Protein has %d residues.\n", protein_last_residue);

        sel2_start = Mol->residue_pointer[protein_last_residue]; 	// first solvent atom
        sel2_end = Mol->N;

        printf("# You selected atoms %d (%s) to %d (%s).\n", sel2_start, Mol->atomnames[sel2_start-1].c_str(), sel2_end, Mol->atomnames[sel2_end-1].c_str());
        if (this->binnary){
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], traj_ans);
        }
        else {
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], gzipped);
        }

        break;

    case 4:
        printf("# Define your ligand residue number: ");
        cin >> sel_res;
        sel_start = Mol->residue_pointer[sel_res-1];
        sel_end = Mol->residue_pointer[sel_res]-1;
        printf("# You selected atoms %d (%s) to %d (%s)\n", sel_start, Mol->atomnames[sel_start-1].c_str(), sel_end, Mol->atomnames[sel_end-1].c_str());
        printf("# Define your second residue: ");
        cin >> sel2_res;
        sel2_start = Mol->residue_pointer[sel2_res-1];
        sel2_end = Mol->residue_pointer[sel2_res]-1;
        printf("# You selected atoms %d (%s) to %d (%s)\n", sel2_start, Mol->atomnames[sel2_start-1].c_str(), sel2_end, Mol->atomnames[sel2_end-1].c_str());
        if (this->binnary){
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], traj_ans);
        }
        else {
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], gzipped);
        }

        break;

    case 5:
        for(unsigned i=0; i<Mol->resnames.size(); i++) {
            if (Mol->resnames[i] != "WAT" and Mol->resnames[i] != "Na+" and Mol->resnames[i] !="Cl-" and Mol->resnames[i] != "Mn2" and Mol->resnames[i] != "Ca2") {
                protein_last_residue++;
            }
        }
        printf("# Protein last residue number: %d.\n", protein_last_residue);

        sel2_start = Mol->residue_pointer[0];			// first solvent atom
        sel2_end = Mol->residue_pointer[protein_last_residue];
        printf("# You selected atoms %d (%s) to %d (%s)\n", sel2_start, Mol->atomnames[sel2_start-1].c_str(), sel2_end, Mol->atomnames[sel2_end-1].c_str());
        if (this->binnary){
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], traj_ans);
        }
        else {
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], gzipped);
        }

        break;

    case 6:
        printf("# Define your first residue range (e.g: 1 100): ");
        cin >> sel_res >> sel2_res ;
        sel_start = Mol->residue_pointer[sel_res-1];
        sel_end = Mol->residue_pointer[sel2_res]-1;
        printf("# You selected atoms %d (%s) to %d (%s)\n", sel_start, Mol->atomnames[sel_start-1].c_str(), sel_end, Mol->atomnames[sel_end-1].c_str());

        printf("# Define your second residue range (e.g: 101 200): ");
        cin >> sel_res >> sel2_res ;
        sel2_start = Mol->residue_pointer[sel_res-1];
        if (sel2_res > Mol->Nres){
            sel2_res = Mol->Nres;
        }
        else if (sel2_res == Mol->Nres){
            sel2_end = Mol->N;
        }
        else {
            sel2_end = Mol->residue_pointer[sel2_res]-1;

            if (sel2_end < 0){
                sel2_end = Mol->N;
            }
        }

        printf("# You selected atoms %d (%s) to %d (%s)\n", sel2_start, Mol->atomnames[sel2_start-1].c_str(), sel2_end, Mol->atomnames[sel2_end-1].c_str());
        if (this->binnary){
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], traj_ans);
        }
        else {
            Coord = new COORD(Mol, sel_start -1, sel_end-1, sel2_start-1, sel2_end-1, argv[2], gzipped);
        }

        break;

    case 7:
        printf("# Define the number of chuncks in the receptor (e.g. 2): ");
        cin >> receptor_chuncks;
        for (int i=0; i< receptor_chuncks; i++){
            printf("# Define the residue range for chunck %2d (e.g: 1 100): ", i);
            cin >> sel_res >> sel2_res ;
            vitmp.push_back(sel_res);
            vitmp.push_back(sel2_res);
            receptor_residues_chuncks.push_back(vitmp);
            vitmp.clear();

            sel_start = Mol->residue_pointer[sel_res-1];
            sel_end = Mol->residue_pointer[sel2_res]-1;
            printf("# You selected atoms %d (%s) to %d (%s)\n", sel_start, Mol->atomnames[sel_start-1].c_str(), sel_end, Mol->atomnames[sel_end-1].c_str());
            vitmp.push_back(sel_start-1);
            vitmp.push_back(sel_end-1);
            receptor_atom_chuncks.push_back(vitmp);
            vitmp.clear();
        }

        printf("# Define the number of chuncks in the ligand (e.g. 2): ");
        cin >> ligand_chuncks;
        for (int i=0; i< ligand_chuncks; i++){
            printf("# Define the residue range for chunck %2d (e.g: 101 200): ", i);
            cin >> sel_res >> sel2_res ;
            vitmp.push_back(sel_res);
            vitmp.push_back(sel2_res);
            ligand_residues_chuncks.push_back(vitmp);
            vitmp.clear();

            sel_start = Mol->residue_pointer[sel_res-1];
            sel_end = Mol->residue_pointer[sel2_res]-1;
            printf("# You selected atoms %d (%s) to %d (%s)\n", sel_start, Mol->atomnames[sel_start-1].c_str(), sel_end, Mol->atomnames[sel_end-1].c_str());
            vitmp.push_back(sel_start-1);
            vitmp.push_back(sel_end-1);
            ligand_atom_chuncks.push_back(vitmp);
            vitmp.clear();
        }

        if (this->binnary){
            Coord = new COORD(Mol, receptor_atom_chuncks, ligand_atom_chuncks, argv[2], traj_ans);
        }
        else {
            printf("Sorry! This option is not available for non binnary trajectories!");
            exit(1);
        }

        break;

    default:
        printf("Selection %d invalid. Please try again.\n", calc_ans);
        break;
    }

    return 0;

}
