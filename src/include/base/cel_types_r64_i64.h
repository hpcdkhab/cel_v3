/*     
 * File:   cel_types.h
 * Project CRESTA (see details on https://cresta-project.eu) Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/
#ifndef CEL_TYPES_H
#define CEL_TYPES_H

#define REAL double
#define INDEX long long
#define COORD double
#define RANK int
#define BOOLEAN int

#define INDEX_MAX 9223372036854775808LL
#define REAL_UNDEF  1.79769e+308
#define EPSILON 1.0e-15


#define INDEX_FORMAT "%6lu"
#define REAL_FORMAT "%g"
#define COORD_FORMAT "%g"
#define RANK_FORMAT "%6d"


#define Nullptr(type) (type*)0
#define POSIX_MEM_ALIGN 32
#endif CEL_TYPES_H
