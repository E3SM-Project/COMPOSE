#ifndef INCLUDE_MPAS_MACROS_HPP
#define INCLUDE_MPAS_MACROS_HPP

#ifdef _OPENMP
# define ompparfor  _Pragma("omp parallel for")
# define omppar     _Pragma("omp parallel")
# define ompfor     _Pragma("omp for")
# define ompbarrier _Pragma("omp barrier")
# define ompsingle  _Pragma("omp single")
# define ompmaster  _Pragma("omp master")
# define ompatomic  _Pragma("omp atomic")
#else
# define ompparfor
# define omppar
# define ompfor
# define ompbarrier
# define ompsingle
# define ompmaster
# define ompatomic
#endif

#endif
