//=============================================================================
//    main.c : implementation
//    Copyright (C) 2018  Bruno Toupance <bruno.toupance@mnhn.fr>
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//=============================================================================

//=============================================================================
#include <stdio.h>
#include <stdlib.h>

#include "d_util.h"


//-----------------------------------------------------------------------------
#include "genpop.h"

//=============================================================================
//  Functions prototypes
//=============================================================================

void  print_usage(char* ExeName);

void  simulate_WrightFisher        (d_size NombreIndividu, d_size NombreGeneration, double* vector_SumHet, double* vector_SumSquareHet, d_size* vector_PopSize);
void  simulate_SexRatio            (d_size NombreMale, d_size NombreFemelle, d_size NombreGeneration, double* vector_SumHet, double* vector_SumSquareHet, d_size* vector_PopSize);
void  simulate_Autofecondation     (d_size NombreIndividu, double TauxAutofecondation, d_size NombreGeneration, double* vector_SumHet, double* vector_SumSquareHet, d_size* vector_PopSize);
void  simulate_CycleDemographique  (d_size NombreIndividu, double TauxCroissance, d_size NombreGenerationParCycle, d_size NombreCycle, double* vector_SumHet, double* vector_SumSquareHet, d_size* vector_PopSize);
void  simulate_Variance            (d_size NombreIndividu, d_size NombreIndividu0, d_size NombreGeneration, double* vector_SumHet, double* vector_SumSquareHet, d_size* vector_PopSize);
//void  simulate_VarianceS           (d_size NombreIndividu, d_size NombreIndividu1, d_size NombreGeneration, double* vector_SumHet, double* vector_SumSquareHet, d_size* vector_PopSize);



//=============================================================================
//  Functions
//=============================================================================


//-----------------------------------------------------------------------------
//  print_usage
//-----------------------------------------------------------------------------
void
print_usage(
	char* ExeName
	){
	(void)fprintf(stderr, "#------------------------------------------------------------------------------\n");
	(void)fprintf(stderr, "# %s: usage\n", ExeName);
	(void)fprintf(stderr, "#------------------------------------------------------------------------------\n");
	(void)fprintf(stderr, "\n");
	(void)fprintf(stderr, "Wright-Fisher:\n");
	(void)fprintf(stderr, "  %s WrightFisher NombreIndividu NombreGeneration NombreSimulation\n", ExeName);
	(void)fprintf(stderr, "\n");
	(void)fprintf(stderr, "Sex-Ratio:\n");
	(void)fprintf(stderr, "  %s SexRatio NombreMale NombreFemelle NombreGeneration NombreSimulation\n", ExeName);
	(void)fprintf(stderr, "\n");
	(void)fprintf(stderr, "Autofecondation:\n");
	(void)fprintf(stderr, "  %s Autofecondation NombreIndividu TauxAutofecondation NombreGeneration NombreSimulation\n", ExeName);
	(void)fprintf(stderr, "\n");
	(void)fprintf(stderr, "Cycle Demographique:\n");
	(void)fprintf(stderr, "  %s CycleDemographique NombreIndividu TauxCroissance NombreGenerationParCycle NombreCycle NombreSimulation\n", ExeName);
	(void)fprintf(stderr, "\n");
	(void)fprintf(stderr, "Variance:\n");
	(void)fprintf(stderr, "  %s Variance NombreIndividu NombreIndividu0 NombreGeneration NombreSimulation\n", ExeName);
	(void)fprintf(stderr, "\n");
}


//=============================================================================
//  Main function
//=============================================================================
int
main(
	int argc,
	char** argv
	){
//-----------------------------------------------------------------------------
	d_size NombreSimulation;
	d_size Simulation;
	d_size NombreGeneration;
	d_size Generation;
	d_size NombreIndividu;
	d_size NombreIndividu0;
	d_size NombreMale;
	d_size NombreFemelle;
	double TauxAutofecondation;
	double TauxCroissance;
	d_size NombreGenerationParCycle;
	d_size NombreCycle;
	T_locus*  loc1;
	d_size loc_idx1;
//-----------------------------------------------------------------------------
	double*   vector_SumHet;
	double*   vector_SumSquareHet;
	d_size*   vector_PopSize;
//-----------------------------------------------------------------------------
	d_byte SimulationModel;
	d_size i;
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
	NombreIndividu = 0;
	NombreIndividu0 = 0;
	NombreMale = 0;
	NombreFemelle = 0;
	NombreGeneration = 0;
	TauxAutofecondation = 0.0;
	TauxCroissance = 1.0;
	NombreGenerationParCycle = 0;
	NombreCycle = 0;
//-----------------------------------------------------------------------------
// Gestion de la ligne de commandes
//-----------------------------------------------------------------------------
	if (argc < 2) {
		print_usage(global__ExeName);
		d_error("Nombre d'arguments incompatible");
	}
//-----------------------------------------------------------------------------
	SimulationModel = 0;
	NombreSimulation = 0;
//-----------------------------------------------------------------------------
//  EXE WrightFisher NombreIndividu NombreGeneration NombreSimulation
//-----------------------------------------------------------------------------
	if (strcmp(argv[1], "WrightFisher") == 0) {
		SimulationModel = 1;
		if (argc != 5) {
			print_usage(global__ExeName);
			d_error("Nombre d'arguments incompatible pour le modele [WrightFisher]");
		}
		i = 2;
		NombreIndividu = atoi(argv[i++]);
		NombreGeneration = atoi(argv[i++]);
		NombreSimulation = atoi(argv[i++]);
	}
//-----------------------------------------------------------------------------
//  EXE SexRatio NombreMale NombreFemelle NombreGeneration NombreSimulation
//-----------------------------------------------------------------------------
	if (strcmp(argv[1], "SexRatio") == 0) {
		SimulationModel = 2;
		if (argc != 6) {
			print_usage(global__ExeName);
			d_error("Nombre d'arguments incompatible pour le modele [SexRatio]");
		}
		i = 2;
		NombreMale = atoi(argv[i++]);
		NombreFemelle = atoi(argv[i++]);
		NombreGeneration = atoi(argv[i++]);
		NombreIndividu = NombreMale+NombreFemelle;
		NombreSimulation = atoi(argv[i++]);
	}
//-----------------------------------------------------------------------------
//  EXE Autofecondation NombreIndividu TauxAutofecondation NombreGeneration NombreSimulation
//-----------------------------------------------------------------------------
	if (strcmp(argv[1], "Autofecondation") == 0) {
		SimulationModel = 3;
		if (argc != 6) {
			print_usage(global__ExeName);
			d_error("Nombre d'arguments incompatible pour le modele [Autofecondation]");
		}
		i = 2;
		NombreIndividu = atoi(argv[i++]);
		TauxAutofecondation = atof(argv[i++]);
		NombreGeneration = atoi(argv[i++]);
		NombreSimulation = atoi(argv[i++]);
	}
//-----------------------------------------------------------------------------
//  EXE CycleDemographique NombreIndividu TauxCroissance NombreGenerationParCycle NombreCycle NombreSimulation
//-----------------------------------------------------------------------------
	if (strcmp(argv[1], "CycleDemographique") == 0) {
		SimulationModel = 4;
		if (argc != 7) {
			print_usage(global__ExeName);
			d_error("Nombre d'arguments incompatible pour le modele [CycleDemographique]");
		}
		i = 2;
		NombreIndividu = atoi(argv[i++]);
		TauxCroissance = atof(argv[i++]);
		NombreGenerationParCycle = atoi(argv[i++]);
		NombreCycle = atoi(argv[i++]);
		NombreSimulation = atoi(argv[i++]);
		NombreGeneration = NombreCycle*NombreGenerationParCycle;
	}
//-----------------------------------------------------------------------------
//  EXE Variance NombreIndividu NombreGeneration NombreSimulation
//-----------------------------------------------------------------------------
	if (strcmp(argv[1], "Variance") == 0) {
		SimulationModel = 5;
		if (argc != 6) {
			print_usage(global__ExeName);
			d_error("Nombre d'arguments incompatible pour le modele [Variance]");
		}
		i = 2;
		NombreIndividu = atoi(argv[i++]);
		NombreIndividu0 = atoi(argv[i++]);
		NombreGeneration = atoi(argv[i++]);
		NombreSimulation = atoi(argv[i++]);
	}
//-----------------------------------------------------------------------------
	if (SimulationModel == 0) {
		print_usage(global__ExeName);
		d_error("Modele inconnu");
	}
//-----------------------------------------------------------------------------
// Random generator initialization
//-----------------------------------------------------------------------------
	global__random = d_random__select("");
//-----------------------------------------------------------------------------
// Genome initialization
//-----------------------------------------------------------------------------
	global__genome = T_genome__new();
//-----------------------------------------------------------------------------
	loc_idx1 = 0;
	loc1 = T_k_alleles__new(loc_idx1, 0.0, 2*NombreIndividu);
	global__genome->m_locus[loc_idx1] = loc1;
//-----------------------------------------------------------------------------
// Heterozygosity and Poplation size
//-----------------------------------------------------------------------------
	vector_SumHet = d_new(double, NombreGeneration+1);
	vector_SumSquareHet = d_new(double, NombreGeneration+1);
	vector_PopSize = d_new(d_size, NombreGeneration+1);
	for (Generation = 0; Generation <= NombreGeneration; Generation++) {
		vector_SumHet[Generation] = 0.0;
		vector_SumSquareHet[Generation] = 0.0;
		vector_PopSize[Generation] = 0;
	}
//-----------------------------------------------------------------------------
// Texte explicatif
//-----------------------------------------------------------------------------
	(void)fprintf(stderr, "#------------------------------------------------------------------------------\n");
	(void)fprintf(stderr, "# Programme [%s]\n", global__ExeName);
	(void)fprintf(stderr, "#------------------------------------------------------------------------------\n");
	(void)fprintf(stderr, "  Ce programme simule la derive genetique.\n");
	(void)fprintf(stderr, "  La simulation est effectuee sur [NombreGeneration] generations non chevauchantes.\n");
	(void)fprintf(stderr, "     NombreGeneration = %ld\n", NombreGeneration);
	(void)fprintf(stderr, "  La simulation est effectuee pour [NombreSimulation] simulations independantes.\n");
	(void)fprintf(stderr, "     NombreSimulation = %ld\n", NombreSimulation);
	(void)fprintf(stderr, "  Les resultats (sortie standard STDOUT) sont les moyennes sur les simulations de l'Heterozygotie a chaque generation.\n");
//-----------------------------------------------------------------------------
	if (SimulationModel == 1) {
		(void)fprintf(stderr, "  Le modele selectionne est : Wright-Fisher\n");
		(void)fprintf(stderr, "  La population est de taille constante [N] individus diploides hermaphrodites.\n");
		(void)fprintf(stderr, "    N = %ld\n", NombreIndividu);
	}
//-----------------------------------------------------------------------------
	if (SimulationModel == 2) {
		(void)fprintf(stderr, "  Le modele selectionne est : Sex-Ratio\n");
		(void)fprintf(stderr, "  La population est de taille constante [N = Nm+Nf] individus diploides males [Nm] ou femelles [Nf].\n");
		(void)fprintf(stderr, "    Nm = %ld\n", NombreMale);
		(void)fprintf(stderr, "    Nf = %ld\n", NombreFemelle);
		(void)fprintf(stderr, "    Sex-Ratio Nm/N = %f\n", NombreMale/((double)(NombreMale+NombreFemelle)));
	}

//-----------------------------------------------------------------------------
	if (SimulationModel == 3) {
		(void)fprintf(stderr, "  Le modele selectionne est : Autofecondation\n");
		(void)fprintf(stderr, "  La population est de taille constante [N] individus diploides hermaphrodites.\n");
		(void)fprintf(stderr, "    N = %ld\n", NombreIndividu);
		(void)fprintf(stderr, "  Une fraction [s] des individus de la population sont produits par autofecondation.\n");
		(void)fprintf(stderr, "    s = %f\n", TauxAutofecondation);
	}

//-----------------------------------------------------------------------------
	if (SimulationModel == 4) {
		(void)fprintf(stderr, "  Le modele selectionne est : CycleDemographique\n");
		(void)fprintf(stderr, "  La population est de taille cyclique [N(t)] individus diploides hermaphrodites.\n");
		(void)fprintf(stderr, "  Au debut de chaque cycle, la taille de la population est N0.\n");
		(void)fprintf(stderr, "    N0 = %ld\n", NombreIndividu);
		(void)fprintf(stderr, "  Le taux de croissance [r] dans un cycle definit N(t+1) = (1+r) N(t).\n");
		(void)fprintf(stderr, "    r = %f\n", TauxCroissance);
		(void)fprintf(stderr, "  Le nombre de cycles est [NombreCycle].\n");
		(void)fprintf(stderr, "    NombreCycle = %ld\n", NombreCycle);
		(void)fprintf(stderr, "  Le nombre de generations par cycle est [NombreGenerationParCycle].\n");
		(void)fprintf(stderr, "    NombreGenerationParCycle = %ld\n", NombreGenerationParCycle);
	}
//-----------------------------------------------------------------------------
	if (SimulationModel == 5) {
		(void)fprintf(stderr, "  Le modele selectionne est : Variance\n");
		(void)fprintf(stderr, "  La population est de taille constante [N] individus diploides hermaphrodites.\n");
		(void)fprintf(stderr, "    N = %ld\n", NombreIndividu);
		(void)fprintf(stderr, "  Le nombre d'individus ne laissant aucun descendant est [N0].\n");
		(void)fprintf(stderr, "    N0 = %ld\n", NombreIndividu0);
		(void)fprintf(stderr, "  La variance du nombre de genes K laisses par chaque individu est :\n");
		(void)fprintf(stderr, "    V[K] = %f\n", 8.0*NombreIndividu0/NombreIndividu);
	}

//-----------------------------------------------------------------------------
// Simulations
//-----------------------------------------------------------------------------
	for (Simulation = 0; Simulation < NombreSimulation; Simulation++) {
//-----------------------------------------------------------------------------
//		(void)fprintf(stderr, "Simulation [%ld]\n", Simulation);
		(void)fprintf(stderr, "."); if ((Simulation+1)%50 == 0) (void)fprintf(stderr, "\n");
//-----------------------------------------------------------------------------
		switch (SimulationModel) {
		case 1:
			simulate_WrightFisher(NombreIndividu, NombreGeneration, vector_SumHet, vector_SumSquareHet, vector_PopSize);
			break;
		case 2:
			simulate_SexRatio(NombreMale, NombreFemelle, NombreGeneration, vector_SumHet, vector_SumSquareHet, vector_PopSize);
			break;
		case 3:
			simulate_Autofecondation(NombreIndividu, TauxAutofecondation, NombreGeneration, vector_SumHet, vector_SumSquareHet, vector_PopSize);
			break;
		case 4:
			simulate_CycleDemographique(NombreIndividu, TauxCroissance, NombreGenerationParCycle, NombreCycle, vector_SumHet, vector_SumSquareHet, vector_PopSize);
			break;
		case 5:
			simulate_Variance(NombreIndividu,NombreIndividu0, NombreGeneration, vector_SumHet, vector_SumSquareHet, vector_PopSize);
			break;
		}
	}
//-----------------------------------------------------------------------------
// Statistiques
//-----------------------------------------------------------------------------
	(void)fprintf(stdout, "%s", "t");
	(void)fprintf(stdout, "\t%s", "PopSize");
	(void)fprintf(stdout, "\t%s", "m_Het");
	(void)fprintf(stdout, "\t%s", "s2_Het");
	(void)fprintf(stdout, "\t%s", "n_Het");
	(void)fprintf(stdout, "\n");
//-----------------------------------------------------------------------------
	for (Generation = 0; Generation <= NombreGeneration; Generation++) {
		(void)fprintf(stdout, "%ld", Generation);
		(void)fprintf(stdout, "\t%ld", vector_PopSize[Generation]);
		(void)fprintf(stdout, "\t%e", vector_SumHet[Generation]/NombreSimulation);
		(void)fprintf(stdout, "\t%e", (vector_SumSquareHet[Generation]-vector_SumHet[Generation]*vector_SumHet[Generation]/NombreSimulation)/(NombreSimulation-1.0));
		(void)fprintf(stdout, "\t%ld", NombreSimulation);
		(void)fprintf(stdout, "\n");
	}
	(void)fprintf(stdout, "\n");
//-----------------------------------------------------------------------------
// Memory managment
//-----------------------------------------------------------------------------
	d_delete(vector_SumHet);
	d_delete(vector_SumSquareHet);
	d_delete(vector_PopSize);
	T_genome__delete(global__genome);
//-----------------------------------------------------------------------------
	return EXIT_SUCCESS;
//-----------------------------------------------------------------------------
}
