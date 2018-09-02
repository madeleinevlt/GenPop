//=============================================================================
//    genepop.c : implementation
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

#define _genpop__implementation_


//=============================================================================
#include <stdio.h>
#include <stdlib.h>

#include "genpop.h"

//-----------------------------------------------------------------------------
// #include ".h"


//=============================================================================
//  T_locus Functions
//=============================================================================

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_locus__delete(
	T_locus*  Locus
){
	d_error_if_fail(Locus != d_null);

	switch (Locus->m_type) {
		case 1: // k alleles
			T_k_alleles__delete(Locus);
			break;
		default:
			d_error("unknown locus type...");
	}
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
double
T_locus__heterozygosity(
	T_locus*       Locus,
	T_population*  Pop
){
	double  Het;

	d_error_if_fail(Locus != d_null);
	d_error_if_fail(Pop != d_null);

	Het = -1.0;

	switch (Locus->m_type) {
		case 1: // k alleles
			Het = T_k_alleles__heterozygosity(Locus, Pop);
			break;
		default:
			d_error("unknown locus type...");
	}

	return Het;
}


//=============================================================================
//  T_k_alleles Functions
//=============================================================================

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_k_alleles__destructor(
	T_k_alleles*  Locus
){
	d_error_if_fail(Locus != d_null);

	d_delete(Locus->m_allele_freq);
	Locus->m_type = 0;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_k_alleles__delete(
	d_ptr  Ptr
){
	T_k_alleles*  Locus;

	d_error_if_fail(Ptr != d_null);

	Locus = (T_k_alleles*)Ptr;

	T_k_alleles__destructor(Locus);

	d_delete(Locus);
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_k_alleles__constructor(
	T_k_alleles*  Locus,
	d_size        idi1,
	double        MutRate,
	d_size        all_cnt1
){
	d_error_if_fail(Locus != d_null);

	Locus->m_type          = 1;  // k alleles
	Locus->m_idi           = idi1;
	Locus->m_mutation_rate = MutRate;

	Locus->m_allele_count = all_cnt1;
	Locus->m_allele_freq = d_new(double, all_cnt1);

	d_memset(Locus->m_allele_freq, 0, Locus->m_allele_count*sizeof(double));
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
T_locus*
T_k_alleles__new(
	d_size    idi1,
	double    MutRate,
	d_size    all_cnt1
){
	T_k_alleles*  Locus;

	Locus = d_new(T_k_alleles, 1);

	T_k_alleles__constructor(Locus, idi1, MutRate, all_cnt1);

	return (T_locus*)Locus;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_k_alleles__frequencies(
	d_ptr          Ptr,
	T_population*  Pop
){
	T_k_alleles*  Locus;
	T_individual* Ind;
	d_size        i;
	d_size        j;

	d_error_if_fail(Ptr != d_null);
	d_error_if_fail(Pop != d_null);
	d_error_if_fail(Pop->m_ind_cnt > 0);

	Locus = (T_k_alleles*)Ptr;

	d_memset(Locus->m_allele_freq, 0, Locus->m_allele_count*sizeof(double));
	for (i = 0; i < Pop->m_ind_cnt; i++) {
		Ind = Pop->m_ind_val+i;
		j = (d_size)(Ind->m_pat.m_genes[Locus->m_idi]);
		Locus->m_allele_freq[j]++;
		j = (d_size)(Ind->m_mat.m_genes[Locus->m_idi]);
		Locus->m_allele_freq[j]++;
	}

	for (j = 0; j < Locus->m_allele_count; j++) {
		Locus->m_allele_freq[j]/= 2.0*Pop->m_ind_cnt;
	}

}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
double
T_k_alleles__heterozygosity(
	d_ptr          Ptr,
	T_population*  Pop
){
	T_k_alleles*  Locus;
	d_size        i;
	double        Het;

	d_error_if_fail(Ptr != d_null);
	d_error_if_fail(Pop != d_null);

	Locus = (T_k_alleles*)Ptr;

	T_k_alleles__frequencies(Locus, Pop);
	Het = 0.0;
	for (i = 0; i < Locus->m_allele_count; i++) {
		Het+= Locus->m_allele_freq[i]*Locus->m_allele_freq[i];
	}
	Het = 1.0-Het;

	return Het;
}



//=============================================================================
//  T_genome Functions
//=============================================================================

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_genome__destructor(
	T_genome*  Genome
){
	d_size i_Locus;

	d_error_if_fail(Genome != d_null);

	for (i_Locus = 0; i_Locus < global__locus_count; i_Locus++) {
		T_locus__delete(Genome->m_locus[i_Locus]);
	}
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_genome__delete(
	T_genome*  Genome
){
	d_error_if_fail(Genome != d_null);

	T_genome__destructor(Genome);

	d_delete(Genome);
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_genome__constructor(
	T_genome*  Genome
){
	d_size     i_Locus;

	d_error_if_fail(Genome != d_null);

	for (i_Locus = 0; i_Locus < global__locus_count; i_Locus++) {
		Genome->m_locus[i_Locus] = d_null;
	}
}


//-----------------------------------------------------------------------------
//  Alloue un nouveau genome
//-----------------------------------------------------------------------------
T_genome*
T_genome__new(
	void
){
	T_genome*  Genome;

	Genome = d_new(T_genome, 1);

	T_genome__constructor(Genome);

	return Genome;
}


//=============================================================================
//  T_gamete Functions
//=============================================================================

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_gamete__destructor(
	T_gamete*  Gamete
){
	d_error_if_fail(Gamete != d_null);
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_gamete__constructor(
	T_gamete*  Gamete
){
	d_error_if_fail(Gamete != d_null);

	d_memset(&(Gamete->m_genes), 0, global__locus_count*sizeof(int));
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_gamete__copy_constructor(
	T_gamete*  dst_Gamete,
	T_gamete*  src_Gamete
){
	d_error_if_fail(dst_Gamete != d_null);
	d_error_if_fail(src_Gamete != d_null);
	d_error_if_fail(src_Gamete != dst_Gamete);

	d_memcpy(&(dst_Gamete->m_genes), &(src_Gamete->m_genes), global__locus_count*sizeof(int));
}


//-----------------------------------------------------------------------------
//  Choisit au hasard le gamete paternel ou le gamete maternel
//-----------------------------------------------------------------------------
void
T_gamete__meiosis_constructor(
	T_gamete*  Gamete,
	T_gamete*  src_Gamete1,
	T_gamete*  src_Gamete2
){
	d_error_if_fail(Gamete != d_null);
	d_error_if_fail(src_Gamete1 != d_null);
	d_error_if_fail(src_Gamete2 != d_null);

	if (d_random__unif(global__random) < 0.5) {
		T_gamete__copy_constructor(Gamete, src_Gamete1);
	} else {
		T_gamete__copy_constructor(Gamete, src_Gamete2);
	}
}


//=============================================================================
//  T_individual Functions
//=============================================================================

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_individual__destructor(
	T_individual*  Ind
){
	d_error_if_fail(Ind != d_null);

	T_gamete__destructor(&(Ind->m_pat));
	T_gamete__destructor(&(Ind->m_mat));

	Ind->m_sex = global__herma;
	Ind->m_gene = 0;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_individual__constructor(
	T_individual*  Ind
){
	d_error_if_fail(Ind != d_null);

	T_gamete__constructor(&(Ind->m_pat));
	T_gamete__constructor(&(Ind->m_mat));
	Ind->m_sex = global__herma;
	Ind->m_gene = 0;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_individual__copy_constructor(
	T_individual*  dst_ind1,
	T_individual*  src_ind1
){

	d_error_if_fail(dst_ind1 != d_null);
	d_error_if_fail(src_ind1 != d_null);
	d_error_if_fail(src_ind1 != dst_ind1);

	T_gamete__copy_constructor(&(dst_ind1->m_pat), &(src_ind1->m_pat));
	T_gamete__copy_constructor(&(dst_ind1->m_mat), &(src_ind1->m_mat));
	dst_ind1->m_sex = src_ind1->m_sex;
	dst_ind1->m_gene = src_ind1->m_gene;

}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_individual__fecondation_constructor(
	T_individual*  Ind,
	T_individual*  Fat,
	T_individual*  Mot,
	d_byte         Sex
){

	d_error_if_fail(Ind != d_null);
	d_error_if_fail(Fat != d_null);
	d_error_if_fail(Mot != d_null);

	T_individual__gamete_get(Fat, &(Ind->m_pat));
	T_individual__gamete_get(Mot, &(Ind->m_mat));
	Ind->m_sex = Sex;
	Ind->m_gene = 0;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_individual__gamete_get(
	T_individual*  Ind,
	T_gamete*      Gamete
){
	d_error_if_fail(Ind != d_null);
	d_error_if_fail(Gamete != d_null);

	T_gamete__meiosis_constructor(Gamete, &(Ind->m_pat), &(Ind->m_mat));
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_individual__print(
	T_individual*  Ind,
	d_size         i_Locus
){
	d_error_if_fail(Ind != d_null);
	d_error_if_fail(i_Locus < global__locus_count);

	(void)fprintf(stderr, "ind: %d %d %X\n", Ind->m_pat.m_genes[i_Locus], Ind->m_mat.m_genes[i_Locus], Ind->m_sex);
}



//=============================================================================
//  T_population Functions
//=============================================================================

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_population__clear(
	T_population*  Pop
){
	d_size  i_Ind;

	d_error_if_fail(Pop != d_null);

	for (i_Ind = 0; i_Ind < Pop->m_ind_cnt; i_Ind++) {
		T_individual__destructor(Pop->m_ind_val+i_Ind);
	}
	Pop->m_ind_cnt = 0;
	Pop->m_mal_cnt = 0;
	Pop->m_fem_cnt = 0;

}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_population__destructor(
	T_population*  Pop
){

	d_error_if_fail(Pop != d_null);

	T_population__clear(Pop);

	d_delete(Pop->m_fem_ptr);
	d_delete(Pop->m_mal_ptr);
	d_delete(Pop->m_ind_val);

	Pop->m_fem_ptr = d_null;
	Pop->m_mal_ptr = d_null;
	Pop->m_ind_val = d_null;

	Pop->m_res_siz = 0;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_population__delete(
	T_population*  Pop
){
	d_error_if_fail(Pop != d_null);

	T_population__destructor(Pop);

	d_delete(Pop);
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_population__constructor(
	T_population*  Pop,
	d_size         Res
){
	d_error_if_fail(Pop != d_null);

// Reserved size
	Pop->m_res_siz = Res;

// Memory allocation
	Pop->m_ind_val = d_new(T_individual, Pop->m_res_siz);
	Pop->m_mal_ptr = d_new(T_individual*, Pop->m_res_siz);
	Pop->m_fem_ptr = d_new(T_individual*, Pop->m_res_siz);

// Actual number of individuals
	Pop->m_ind_cnt = 0;
	Pop->m_fem_cnt = 0;
	Pop->m_mal_cnt = 0;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
T_population*
T_population__new(
	d_size  Res
){
	T_population*  Pop;

	Pop = d_new(T_population, 1);

	T_population__constructor(Pop, Res);

	return Pop;
}


//-----------------------------------------------------------------------------
//  Ajoute une copie de l'individu [Ind] dans la population [Pop]
//-----------------------------------------------------------------------------
T_individual*
T_population__ind_add(
	T_population*  Pop,
	T_individual*  Ind
){
	d_size  i_Ind;

	d_error_if_fail(Pop != d_null);  // Verification de la validite de la [Pop]
	d_error_if_fail(Pop->m_res_siz >= Pop->m_ind_cnt);  // Verification de la disponibilite d'espace alloue

	i_Ind = Pop->m_ind_cnt;  // Indice du nouvel individu dans le tableau des individus
	T_individual__copy_constructor(Pop->m_ind_val+i_Ind, Ind);  // Copie des donnees de l'individu
	Pop->m_ind_cnt++;  // Incrementation du nombre d'individus

	if (Ind->m_sex & global__male) {  // Si l'individu possede la fonction male
		Pop->m_mal_ptr[Pop->m_mal_cnt] = Pop->m_ind_val+i_Ind;  // Ajout d'un pointeur
		Pop->m_mal_cnt++;  // Incrementation du nombre de males
	}
	if (Ind->m_sex & global__female) {  // Si l'individu possede la fonction femelle
		Pop->m_fem_ptr[Pop->m_fem_cnt] = Pop->m_ind_val+i_Ind;  // Ajout d'un pointeur
		Pop->m_fem_cnt++;  // Incrementation du nombre de femelles
	}
	return(Pop->m_ind_val+i_Ind);
}


//-----------------------------------------------------------------------------
// Cree un individu de sexe [Sex] a partir du pere [Fat] et de la mere [Mot]
// Ajoute cet individu a la population [Pop]
//-----------------------------------------------------------------------------
T_individual*
T_population__child_add(
	T_population*  Pop,
	T_individual*  Fat,
	T_individual*  Mot,
	d_byte         Sex
){
	T_individual  Ind;

	d_error_if_fail(Pop != d_null);  // Verification de la validite de la [Pop]
	d_error_if_fail(Fat != d_null);  // Verification de la validite du pere [Fat]
	d_error_if_fail(Mot != d_null);  // Verification de la validite de la mere [Mot]

	T_individual__fecondation_constructor(&Ind, Fat, Mot, Sex);  // Creation de l'individu
	return(T_population__ind_add(Pop, &Ind));  // Ajout de l'individu
}


//-----------------------------------------------------------------------------
// Tire un individu de sexe [Sex] dans la population [Pop]
// Tous les individus de sexe [Sex] de la population [Pop] sont equiprobables
// Renvoie un pointeur sur cet individu
//-----------------------------------------------------------------------------
T_individual*
T_population__ind_get(
	T_population*  Pop,
	d_byte         Sex
){
	d_size         idx1;
	T_individual*  Ind;

	d_error_if_fail(Pop != d_null);

	switch (Sex) {
		case global__male:
			if (Pop->m_mal_cnt == 0) {
				d_error("no males...");
			}
			idx1 = (d_size)(d_random__unif(global__random)*Pop->m_mal_cnt);
			Ind = Pop->m_mal_ptr[idx1];
			break;

		case global__female:
			if (Pop->m_fem_cnt == 0) {
				d_error("no females...");
			}
			idx1 = (d_size)(d_random__unif(global__random)*Pop->m_fem_cnt);
			Ind = Pop->m_fem_ptr[idx1];
			break;

		default:
			if (Pop->m_ind_cnt == 0) {
				d_error("no individuals...");
			}
			idx1 = (d_size)(d_random__unif(global__random)*Pop->m_ind_cnt);
			Ind = Pop->m_ind_val+idx1;
	}

	if (Ind == d_null) {
		d_error("null individual...");
	}

	return Ind;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
d_boolean
T_population__is_fixed(
	T_population*  Pop,
	d_size         i_Locus,
	int*           all_sta1
){
	d_size         i_Ind;
	d_boolean      fix1;
	T_individual*  Ind;

	d_error_if_fail(Pop != d_null);
	d_error_if_fail(i_Locus < global__locus_count);
	d_error_if_fail(Pop->m_ind_cnt > 0);

	Ind = Pop->m_ind_val;
	*all_sta1 = Ind->m_pat.m_genes[i_Locus];

	i_Ind = 0;
	fix1 = d_true;

	while ( fix1 && (i_Ind < Pop->m_ind_cnt) ) {
		Ind = Pop->m_ind_val+i_Ind;
		fix1 = (*all_sta1 == Ind->m_pat.m_genes[i_Locus]) && (*all_sta1 == Ind->m_mat.m_genes[i_Locus]);
		i_Ind++;
	}

	return fix1;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
T_population__print(
	T_population*  Pop,
	d_size         i_Locus
){
	d_size         i_Ind;
	T_individual*  Ind;

	d_error_if_fail(Pop != d_null);
	d_error_if_fail(i_Locus < global__locus_count);
	d_error_if_fail(Pop->m_ind_cnt > 0);


	(void)fprintf(stderr, "[%ld][%ld:%ld:%ld]", Pop->m_res_siz, Pop->m_ind_cnt, Pop->m_mal_cnt, Pop->m_fem_cnt);
	for (i_Ind = 0; i_Ind < Pop->m_ind_cnt; i_Ind++) {
		Ind = Pop->m_ind_val+i_Ind;
		(void)fprintf(stderr, " %02d %02d", Ind->m_pat.m_genes[i_Locus], Ind->m_mat.m_genes[i_Locus]);
	}
	(void)fprintf(stderr, " %f", T_population__heterozygosity(Pop, i_Locus));
	(void)fprintf(stderr, "\n");
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
double
T_population__heterozygosity(
	T_population*  Pop,
	d_size         i_Locus
){

	d_error_if_fail(Pop != d_null);
	d_error_if_fail(i_Locus < global__locus_count);
	d_error_if_fail(Pop->m_ind_cnt > 0);


	return T_locus__heterozygosity(global__genome->m_locus[i_Locus], Pop);
}



