//=============================================================================
//    genepop.h : interface
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

#ifndef _genpop__interface_
#define _genpop__interface_

//-----------------------------------------------------------------------------
#ifdef _genpop__implementation_
#define _where_genpop_
#else
#define _where_genpop_ extern
#endif  // _genpop__implementation_


//=============================================================================
#include "d_util.h"


//=============================================================================
#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus


//=============================================================================
//  Defines
//=============================================================================
#define global__male    0x01  // binary 01 == male
#define global__female  0x02  // binary 10 == female
#define global__herma   0x03  // binary 11 == hermaphrodite
//-----------------------------------------------------------------------------
#define global__locus_count  1  // number of loci
//-----------------------------------------------------------------------------
#define global__ExeName "mbi_ne"
//-----------------------------------------------------------------------------


//=============================================================================
//  Types
//=============================================================================

typedef struct _T_locus_         T_locus;
typedef struct _T_k_alleles_     T_k_alleles;

typedef struct _T_genome_        T_genome;

typedef struct _T_gamete_        T_gamete;
typedef struct _T_individual_    T_individual;
typedef struct _T_population_    T_population;


//=============================================================================
//  Structures
//=============================================================================

//=============================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_locus
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

struct _T_locus_
{
	int     m_type;           // locus type:
	                          //    1 == k alleles
	d_size  m_idi;            // locus identifier
	double  m_mutation_rate;  // mutation rate
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_k_alleles
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

struct _T_k_alleles_
{
	int      m_type;           // locus type == 1
	d_size   m_idi;            // locus identifier
	double   m_mutation_rate;  // mutation rate

	d_size   m_allele_count;   // number of alleles
	double*  m_allele_freq;    // allele frequencies
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_genome
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

struct _T_genome_
{
	T_locus*  m_locus[global__locus_count];  // loci
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_gamete
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

struct _T_gamete_
{
	int  m_genes[global__locus_count];  // genes
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_individu
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

struct _T_individual_
{
	T_gamete  m_pat;    // paternal gamete
	T_gamete  m_mat;    // maternal gamete
	d_byte    m_sex;    // sex of individual
	
	d_size    m_gene;   // number of genes
	int       m_var;    // number of genes to transmit
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_population
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

struct _T_population_
{
	d_size          m_res_siz;  // reserved size to store individuals

	d_size          m_ind_cnt;  // actual number of individuals
	T_individual*   m_ind_val;  // array of individuals

	d_size          m_mal_cnt;  // actual number of males
	T_individual**  m_mal_ptr;  // array of pointers on males

	d_size          m_fem_cnt;  // actual number of females
	T_individual**  m_fem_ptr;  // array of pointers on females
};


//=============================================================================
//  Variables
//=============================================================================

_where_genpop_  T_genome*  global__genome;
_where_genpop_  d_random*  global__random;


//=============================================================================
//  Functions prototypes
//=============================================================================

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_locus
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_where_genpop_ void       T_locus__delete         (T_locus* Locus);

_where_genpop_ double     T_locus__heterozygosity (T_locus* Locus, T_population* Pop);


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_k_alleles
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_where_genpop_ void          T_k_alleles__destructor     (T_k_alleles* Locus);
_where_genpop_ void          T_k_alleles__delete         (d_ptr Ptr);

_where_genpop_ void          T_k_alleles__constructor    (T_k_alleles* Locus, d_size idi1, double MutRate, d_size all_cnt1);
_where_genpop_ T_locus*      T_k_alleles__new            (d_size idi1, double MutRate, d_size all_cnt1);

_where_genpop_ void          T_k_alleles__frequencies    (d_ptr Ptr, T_population* Pop);
_where_genpop_ double        T_k_alleles__heterozygosity (d_ptr Ptr, T_population* Pop);


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_genome
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_where_genpop_ void       T_genome__destructor  (T_genome* Genome);
_where_genpop_ void       T_genome__delete      (T_genome* Genome);

_where_genpop_ void       T_genome__constructor (T_genome* Genome);
_where_genpop_ T_genome*  T_genome__new         (void);

_where_genpop_ void       T_genome__locus_push  (T_genome* Genome, T_locus* Locus);
_where_genpop_ d_size     T_genome__locus_count (T_genome* Genome);


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_gamete
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_where_genpop_ void       T_gamete__destructor          (T_gamete* Gamete);

_where_genpop_ void       T_gamete__constructor         (T_gamete* Gamete);
_where_genpop_ void       T_gamete__copy_constructor    (T_gamete* dst_Gamete, T_gamete* src_Gamete);
_where_genpop_ void       T_gamete__meoisis_constructor (T_gamete* Gamete, T_gamete* src_Gamete1, T_gamete* src_Gamete2);


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_individual
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_where_genpop_ void           T_individual__destructor              (T_individual* Ind);

_where_genpop_ void           T_individual__constructor             (T_individual* Ind);
_where_genpop_ void           T_individual__copy_constructor        (T_individual* dst_Ind, T_individual* src_Ind);
_where_genpop_ void           T_individual__fecondation_constructor (T_individual* Ind, T_individual* Fat, T_individual* Mot, d_byte Sex);

_where_genpop_ void           T_individual__gamete_get              (T_individual* Ind, T_gamete* Gamete);

_where_genpop_ void           T_individual__print                   (T_individual* Ind, d_size i_Locus);


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  T_population
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_where_genpop_ void           T_population__destructor     (T_population* Pop);
_where_genpop_ void           T_population__delete         (T_population* Pop);
_where_genpop_ void           T_population__clear          (T_population* Pop);

_where_genpop_ void           T_population__constructor    (T_population* Pop, d_size Res);
_where_genpop_ T_population*  T_population__new            (d_size Res);

_where_genpop_ T_individual*  T_population__ind_add        (T_population* Pop, T_individual* Ind);
_where_genpop_ T_individual*  T_population__child_add      (T_population* Pop, T_individual* Fat, T_individual* Mot, d_byte Sex);

_where_genpop_ T_individual*  T_population__ind_get        (T_population* Pop, d_byte Sex);

_where_genpop_ d_boolean      T_population__is_fixed       (T_population* Pop, d_size i_Locus, int* all_sta1);
_where_genpop_ void           T_population__print          (T_population* Pop, d_size i_Locus);

_where_genpop_ double         T_population__heterozygosity (T_population* Pop, d_size i_Locus);


#ifdef __cplusplus
}
#endif  // __cplusplus

#endif  // _genpop__interface_


