//=============================================================================
//    d_util.h : interface
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

#ifndef _d_util__interface_
#define _d_util__interface_

//-----------------------------------------------------------------------------
#ifdef _d_util__implementation_
#define _where__d_util_
#else
#define _where__d_util_ extern
#endif  // _d_util__implementation_


//=============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


//=============================================================================
#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus



//=============================================================================
//  Constants
//=============================================================================

#ifndef d_E
#define d_E        2.71828182845904523536028747135  // e
#endif

#ifndef d_PI
#define d_PI       3.14159265358979323846264338328  // pi
#endif

#ifndef d_SQRT2
#define d_SQRT2    1.41421356237309504880168872421  // sqrt(2)
#endif

#ifndef d_EULER
#define d_EULER    0.57721566490153286060651209008  // Euler constant
#endif


//=============================================================================
//  Defines & Types
//=============================================================================

//-----------------------------------------------------------------------------
//  Generic types
//-----------------------------------------------------------------------------

#ifdef	NULL  // null pointer
#define	d_null  NULL
#else
#define	d_null  ((void*) 0)
#endif

#define d_max(a1, a2)  (((a1) > (a2)) ?  (a1) : (a2))  // maximum
#define d_min(a1, a2)  (((a1) < (a2)) ?  (a1) : (a2))  // minimum
#define d_abs(a1)      (((a1) < 0  ) ? -(a1) : (a1))   // absolute value


typedef unsigned char   d_byte;     // byte
typedef char            d_char;     // character
typedef int             d_int;      // integer

typedef d_int           d_boolean;  // boolean
#define d_false         0           // false
#define d_true          1           // true

typedef void*           d_ptr;      // generic pointer

typedef unsigned long   d_size;     // interger for memory access


//-----------------------------------------------------------------------------
//  Random types
//-----------------------------------------------------------------------------
typedef struct _d_random_  d_random;

struct _d_random_ {
	void      (*m_set)(void* RandomGeneratorState, long int RandomGeneratorSeed);
	long int  (*m_get)(void* RandomGeneratorState);
	double    (*m_get_double)(void* RandomGeneratorState);
	void*       m_state;
};


//=============================================================================
//  Memory management functions prototypes
//=============================================================================

_where__d_util_  d_ptr  d_malloc           (d_size siz1);
_where__d_util_  d_ptr  d_realloc          (d_ptr ptr1, d_size siz1);
_where__d_util_  void   d_free             (d_ptr ptr1);
_where__d_util_  d_ptr  d_memmove          (d_ptr dst1, d_ptr src1, d_size siz1);
_where__d_util_  d_ptr  d_memcpy           (d_ptr dst1, d_ptr src1, d_size siz1);
_where__d_util_  d_ptr  d_memset           (d_ptr dst1, d_int val1, d_size siz1);


#define d_new(T_val, val_cnt1) \
	((T_val*)d_malloc((d_size)sizeof(T_val)*(val_cnt1)))

#define d_resize(T_val, ptr1, val_cnt1) \
	((T_val*)d_realloc(ptr1, (d_size)sizeof(T_val)*(val_cnt1)))

#define d_delete \
	d_free


//=============================================================================
//  Error functions prototypes
//=============================================================================

_where__d_util_  void  d_error              (d_char*  ErrorMessage);

#define d_error_if_fail(TestExpression) \
	if (!(TestExpression)) { \
		(void)fprintf(stderr, \
		              "\nERROR: assertion '%s' failed.\n", \
		              #TestExpression); \
		exit(EXIT_FAILURE); \
	};


//=============================================================================
//  d_random function prototypes
//=============================================================================

_where__d_util_  d_random*  d_random__select       (char* RandomGeneratorName);
_where__d_util_  void       d_random__seed         (d_random* RandomGenerator, long int RandomGeneratorSeed);
_where__d_util_  void       d_random__default_seed (d_random* RandomGenerator);

_where__d_util_  double     d_random__unif         (d_random* RandomGenerator);
_where__d_util_  double     d_random__norm         (d_random* RandomGenerator);
_where__d_util_  double     d_random__expo         (d_random* RandomGenerator);

_where__d_util_  void       d_random__permutation         (d_random* RandomGenerator, d_size* Vector, d_size Length);
_where__d_util_  void       d_random__permutation_d_byte  (d_random* RandomGenerator, d_byte* Vector, d_size Length);
_where__d_util_  void       d_random__permutation_int     (d_random* RandomGenerator, int* Vector, d_size Length);


#ifdef __cplusplus
}
#endif  // __cplusplus

#endif  // _d_util__interface_

