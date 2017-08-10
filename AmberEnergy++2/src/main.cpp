/*
 * main.cpp
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

#include<stdlib.h>
#include<stdio.h>
#include"RunEngine.h"

using namespace std;

int main (int argc, char* argv[]) {

	if (argc < 2) {
                printf("USAGE: %s prmtop mdcrd\n", argv[0]);
		exit(1);
	}

	RunEngine Engine;
	Engine.Run(argv);

}
