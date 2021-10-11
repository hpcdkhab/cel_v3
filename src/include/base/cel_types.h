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
#include <limits.h>
#define REAL double
#define INDEX int//long long//int
#define COORD double
#define RANK int
#define BOOLEAN int

#define INDEX_MAX INT_MAX //9223372036854775807LL//4294967295//LONG_MAX
#define REAL_UNDEF  1.79769e+308
#define EPSILON 1.0e-15


#define INDEX_FORMAT "%6lld"
#define REAL_FORMAT "%g"
#define COORD_FORMAT "%g"
#define RANK_FORMAT "%6d"


#define Nullptr(type) (type*)0
#define POSIX_MEM_ALIGN 32
#endif //CEL_TYPES_H
