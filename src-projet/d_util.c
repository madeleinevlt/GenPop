//=============================================================================
//    d_util.c : implementation
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

#define _d_util__implementation_


//=============================================================================
#include <time.h>

//-----------------------------------------------------------------------------
#include "d_util.h"


//=============================================================================
//  Memory managment functions
//=============================================================================
#undef D_MEM_TRACE

//-----------------------------------------------------------------------------
//  d_malloc
//-----------------------------------------------------------------------------
d_ptr
d_malloc(
	d_size  siz1
){
	d_ptr  ptr0;

	if (siz1 == 0) {
		return d_null;
	}

	ptr0 = (d_ptr)malloc(siz1);
	if (ptr0 == d_null) {
		d_error("could not allocate");
	}

#ifdef D_MEM_TRACE
	(void)fprintf(stderr, "%p\t+1\n", ptr0);
#endif

	return ptr0;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
d_ptr
d_realloc(
	d_ptr   ptr1,
	d_size  siz1
){
	d_ptr  ptr0;

	if (siz1 == 0) {
		d_free(ptr1);

#ifdef D_MEM_TRACE
	(void)fprintf(stderr, "%p\t-1\n", ptr1);
#endif

		return d_null;
	}

	if (ptr1 == d_null) {
		ptr0 = (d_ptr)malloc(siz1);

#ifdef D_MEM_TRACE
	(void)fprintf(stderr, "%p\t+1\n", ptr0);
#endif

	} else {
		ptr0 = (d_ptr)realloc(ptr1, siz1);

#ifdef D_MEM_TRACE
	(void)fprintf(stderr, "%p\t-1\n", ptr1);
	(void)fprintf(stderr, "%p\t+1\n", ptr0);
#endif

	}

	if (ptr0 == d_null) {
		d_error("could not reallocate");
	}

	return ptr0;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
d_free(
	d_ptr  ptr1
){
	if (ptr1 != d_null) {
		(void)free(ptr1);

#ifdef D_MEM_TRACE
	(void)fprintf(stderr, "%p\t-1\n", ptr1);
#endif

	}
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
d_ptr
d_memmove(
	d_ptr   dst1,
	d_ptr   src1,
	d_size  siz1
){
	return memmove(dst1, src1, siz1);
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
d_ptr
d_memcpy(
	d_ptr   dst1,
	d_ptr   src1,
	d_size  siz1
){
	return memcpy(dst1, src1, siz1);
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
d_ptr
d_memset(
	d_ptr   dst1,
	d_int   val1,
	d_size  siz1
){
	return memset(dst1, val1, siz1);
}



//=============================================================================
//  Error functions
//=============================================================================

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
d_error(
	d_char*  ErrorMessage
){
	(void)fprintf(stderr, "ERROR: %s\n", ErrorMessage);
	(void)fflush(stderr);
	exit(EXIT_FAILURE);
}




//=============================================================================
//  Random Generator [ran3]
//=============================================================================

//-----------------------------------------------------------------------------
//  Random state structure & type for [ran3]
//-----------------------------------------------------------------------------
typedef struct _d_random__ran3_state_  d_random__ran3_state;

struct _d_random__ran3_state_ {
	unsigned int  m_i;
	unsigned int  m_j;
	long int      m_v[56];
};

//-----------------------------------------------------------------------------
//  Constants for [ran3]
//-----------------------------------------------------------------------------
#define d_random__ran3_BIG  1000000000
#define d_random__ran3_SEED  161803398


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
long int
d_random__ran3_get (
	void* RandomGeneratorState
){
	d_random__ran3_state*  s1 = (d_random__ran3_state*)(RandomGeneratorState);
	long int               v1;

	s1->m_i++;
	if (s1->m_i == 56) {
		s1->m_i = 1;
	}

	s1->m_j++;
	if (s1->m_j == 56) {
		s1->m_j = 1;
	}


	v1 = s1->m_v[s1->m_i]-s1->m_v[s1->m_j];

	if (v1 < 0) {
		v1+= d_random__ran3_BIG;
	}

	s1->m_v[s1->m_i] = v1;

	return v1;
}



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
double
d_random__ran3_get_double(
	void* RandomGeneratorState
){
	return d_random__ran3_get(RandomGeneratorState)/(double)d_random__ran3_BIG ;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void
d_random__ran3_set(
	void*     RandomGeneratorState,
	long int  see1
){
	d_random__ran3_state*  s1 = (d_random__ran3_state*)(RandomGeneratorState);
	unsigned int   i;
	unsigned int   j;
	long int       v1;
	long int       v2;

	if (see1 == 0) {
		see1 = 1;
	}

	v1 = d_abs(d_random__ran3_SEED-d_abs(see1));
	v1%= d_random__ran3_BIG;

	s1->m_v[0] = 0;
	s1->m_v[55] = v1;

	v2 = 1;
	for (i = 1; i < 55; i++) {
		j = (21*i)%55;
		s1->m_v[j] = v2;
		v2 = v1-v2;
		if (v2 < 0) {
			v2+= d_random__ran3_BIG;
		}
		v1 = s1->m_v[j];
	}

	for (j = 0; j < 4; j++) {
		for (i = 1; i < 56; i++) {
			v1 = s1->m_v[i]-s1->m_v[1+(i+30)%55];
			if (v1 < 0) {
				v1+= d_random__ran3_BIG;
			}
			s1->m_v[i] = v1;
		}
	}

	s1->m_i = 0;
	s1->m_j = 31;

	return;
}

//-----------------------------------------------------------------------------
//  State for [ran3]
//-----------------------------------------------------------------------------
static d_random__ran3_state  static__ran3_state;

//-----------------------------------------------------------------------------
//  Random generator for [ran3]
//-----------------------------------------------------------------------------
static d_random  static__ran3 =
{
	d_random__ran3_set,
	d_random__ran3_get,
	d_random__ran3_get_double,
	&static__ran3_state
};




//=============================================================================
//  Generic random funtions for a random generator
//=============================================================================

//-----------------------------------------------------------------------------
//  Returns the random generator [RandomGeneratorName]
//    N.B.: actually, only selects [ran3]
//-----------------------------------------------------------------------------
d_random*
d_random__select(
	char* RandomGeneratorName
){
	d_random*  RandomGenerator;

	RandomGenerator = &static__ran3;
	d_random__default_seed(RandomGenerator);

	return RandomGenerator;
}


//-----------------------------------------------------------------------------
//  Sets the seed for the random generator [RandomGenerator]
//-----------------------------------------------------------------------------
void
d_random__seed(
	d_random*  RandomGenerator,
	long int   RandomGeneratorSeed
){
	RandomGenerator->m_set(RandomGenerator->m_state, RandomGeneratorSeed);
}


//-----------------------------------------------------------------------------
//  Sets a default seed for the random generator [RandomGenerator]
//-----------------------------------------------------------------------------
void
d_random__default_seed(
	d_random*  RandomGenerator
){
	time_t  CurrentTime;

	(void)time(&CurrentTime);
	RandomGenerator->m_set(RandomGenerator->m_state, CurrentTime);
}


//-----------------------------------------------------------------------------
//  Return an Uniform deviate on ]0, 1[
//-----------------------------------------------------------------------------
double
d_random__unif(
	d_random* RandomGenerator
){
	double x;

	do {
		x = RandomGenerator->m_get_double(RandomGenerator->m_state);
	} while ((x >= 1.0) || (x <= 0.0));

	return x;
}


//-----------------------------------------------------------------------------
//  Returns a Normal (mu=0, sigma=1) deviate
//-----------------------------------------------------------------------------
double
d_random__norm(
	d_random* RandomGenerator
){
	double  x;
	double  y;
	double  r;
	double  z;

	do {
		x = -1.0+2.0*d_random__unif(RandomGenerator);
		y = -1.0+2.0*d_random__unif(RandomGenerator);
		r = x*x+y*y;
	} while ( (r > 1.0) || (r == 0.0) );

	z = y*sqrt(-2.0*log(r)/r);

	return z;
}


//-----------------------------------------------------------------------------
//  Returns an Exponential (lambda=1) deviate
//-----------------------------------------------------------------------------
double
d_random__expo(
	d_random* RandomGenerator
){
	return -log(d_random__unif(RandomGenerator));
}


//-----------------------------------------------------------------------------
//  Permutation of a d_size vector
//-----------------------------------------------------------------------------
void
d_random__permutation(
	d_random* RandomGenerator,
    d_size*   Vector,
    d_size    Length
){
    d_size Pos;
    d_size i;
    d_size Temp;

    for (i = 0; i < Length-1; i++) {
        Pos = i+(d_size)((Length-i)*d_random__unif(RandomGenerator));
        Temp = Vector[i];
        Vector[i] = Vector[Pos];
        Vector[Pos] = Temp;
    }
}

//-----------------------------------------------------------------------------
//  Permutation of a d_byte vector
//-----------------------------------------------------------------------------
void
d_random__permutation_d_byte(
	d_random* RandomGenerator,
    d_byte*   Vector,
    d_size    Length
){
    d_size Pos;
    d_size i;
    d_byte Temp;

    for (i = 0; i < Length-1; i++) {
        Pos = i+(d_size)((Length-i)*d_random__unif(RandomGenerator));
        Temp = Vector[i];
        Vector[i] = Vector[Pos];
        Vector[Pos] = Temp;
    }
}

//-----------------------------------------------------------------------------
//  Permutation of a int vector
//-----------------------------------------------------------------------------
void
d_random__permutation_int(
	d_random* RandomGenerator,
    int*      Vector,
    d_size    Length
){
    d_size Pos;
    d_size i;
    int    Temp;

    for (i = 0; i < Length-1; i++) {
        Pos = i+(d_size)((Length-i)*d_random__unif(RandomGenerator));
        Temp = Vector[i];
        Vector[i] = Vector[Pos];
        Vector[Pos] = Temp;
    }
}

