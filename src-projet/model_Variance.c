//=============================================================================
//    model_Variance.c : implementation
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
//  Functions
//=============================================================================

//-----------------------------------------------------------------------------
//  simulate_Variance
//-----------------------------------------------------------------------------
void
simulate_Variance(
	d_size NombreIndividu,
	d_size NombreIndividu0,
	d_size NombreGeneration,
	double*  vector_SumHet,
	double*  vector_SumSquareHet,
	d_size*  vector_PopSize
	){

	//-----------------------------------------------------------------------------
	d_size i;
	T_individual Ind;
	d_byte Sex;
	//-----------------------------------------------------------------------------
	d_size Generation;
	//-----------------------------------------------------------------------------
	T_individual*  Father;                  // father
	T_individual*  Mother;                  // mother
	//-----------------------------------------------------------------------------
	T_population*  ParentPop;                 // parent population
	T_population*  OffspringPop;                 // offspring population
	T_population*  TempPop;                 // population pointer
	//-----------------------------------------------------------------------------
	double Het;                          // heterozygosity
	//-----------------------------------------------------------------------------
	// Memory allocation
	//-----------------------------------------------------------------------------
	ParentPop    = T_population__new(NombreIndividu);
	OffspringPop = T_population__new(NombreIndividu);
	//-----------------------------------------------------------------------------
	// Population initialization
	//-----------------------------------------------------------------------------
	for (i = 0; i < NombreIndividu; i++) {
		T_individual__constructor(&Ind);
		//-----------------------------------------------------------------------------
		Ind.m_pat.m_genes[0] = i*2;
		Ind.m_mat.m_genes[0] = i*2+1;
		Ind.m_sex = global__herma;
		//-----------------------------------------------------------------------------
		T_population__ind_add(ParentPop, &Ind);
		//-----------------------------------------------------------------------------
		T_individual__destructor(&Ind);

		//-----------------------------------------------------------------------------
	}

	//-----------------------------------------------------------------------------
	// Simulation
	//-----------------------------------------------------------------------------
	for (Generation = 0; Generation <= NombreGeneration; Generation++) {
		//-----------------------------------------------------------------------------

		Het = T_population__heterozygosity(ParentPop, 0);

		vector_SumHet[Generation]+= Het;
		vector_SumSquareHet[Generation]+= Het*Het;
		//-----------------------------------------------------------------------------
		vector_PopSize[Generation] = NombreIndividu;
		d_size N2 = NombreIndividu - (2*NombreIndividu0);
		d_size N4 = NombreIndividu0;



		//-----------------------------------------------------------------------------
		for (i = 0; i < NombreIndividu; i++) {
			Father = T_population__ind_get(ParentPop, global__male);
			Mother = T_population__ind_get(ParentPop, global__female);
			//-----------------------------------------------------------------------------
			Sex = global__herma;
			//-----------------------------------------------------------------------------
			if (i<=N2) {
				T_population__child_add(OffspringPop, Father, Mother, Sex);
			}
			else if (i>N2 && i<=N4) {
				T_population__child_add(OffspringPop, Father, Mother, Sex);
				Mother = T_population__ind_get(ParentPop, global__female);
				T_population__child_add(OffspringPop, Father, Mother, Sex);

			}
			//-----------------------------------------------------------------------------
		}
		//-----------------------------------------------------------------------------
		T_population__clear(ParentPop);
		//-----------------------------------------------------------------------------
		TempPop      = ParentPop;
		ParentPop    = OffspringPop;
		OffspringPop = TempPop;
		//-----------------------------------------------------------------------------
	}
	//-----------------------------------------------------------------------------
	// Memory managment
	//-----------------------------------------------------------------------------
	T_population__delete(ParentPop);
	T_population__delete(OffspringPop);
	//-----------------------------------------------------------------------------

}
