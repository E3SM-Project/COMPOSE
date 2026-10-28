#ifndef INCLUDE_MPAS_TRIQUAD_HPP
#define INCLUDE_MPAS_TRIQUAD_HPP

#include "mpas_types.hpp"

struct TriangleQuadrature {
  static void get_coef(const int order, const Real*& coord, const Real*& weight, int& n);

private:
  static const Real trisym_order2_coord_  [  9];
  static const Real trisym_order2_weight_ [  3];
  static const Real trisym_order4_coord_  [ 18];
  static const Real trisym_order4_weight_ [  6];
  static const Real tritay_order6_coord_  [ 33];
  static const Real tritay_order6_weight_ [ 11];
  static const Real trisym_order8_coord_  [ 48];
  static const Real trisym_order8_weight_ [ 16];
  static const Real tritay_order12_coord_ [ 96];
  static const Real tritay_order12_weight_[ 32];
  static const Real trisym_order20_coord_ [264];
  static const Real trisym_order20_weight_[ 88];
};

#endif
