#include "mpas_check.hpp"
#include "mpas_geometry.hpp"

void check (const MpasMesh& m) {
  const auto& t = m.t; const auto& c = m.c; const auto& e = m.e;
  const auto& v = m.v; const auto& g = m.g;

  static const Real Rearth = 6.371e6;
  // Normalize all points now that we know they are all on the sphere.
  static const bool nmlz = true;
  Real R; {
    const int ci = 0;
    const Real cc[3] = {c.x[ci], c.y[ci], c.z[ci]};
    if (geo::norm(cc) < 2) R = 1; else R = Rearth;
  }
  
  int out = 10;
  Stats c_rad_s, e_rad_s, v_rad_s, ec_from_v_s, ec_from_c_s, ortho_s,
    dvEdge_s, dcEdge_s, area_s, angle_s, centroid_s;

  for (int ci = 0; ci < t.ncell; ++ci) {
    Vec cc;
    fill(c, ci, cc);
    c_rad_s(std::abs(geo::norm(cc) - R)/R); // Cell center points on the sphere.
    geo::normalize(cc);

    const int ne = t.c2ne[ci];
    for (int i = 0; i < ne; ++i) {
      const int ei = t.get_c2e(ci,i), vi0 = t.get_e2v(ei,0), vi1 = t.get_e2v(ei,1),
        ci0 = t.get_e2c(ei,0), ci1 = t.get_e2c(ei,1);

      { // Edge center points on the sphere.
        Vec ec; fill(e, ei, ec, false);
        e_rad_s(std::abs(geo::norm(ec) - R)/R);
      }

      Vec ec; fill(e, ei, ec, nmlz);
      Vec v0; fill(v, vi0, v0, nmlz);
      Vec v1; fill(v, vi1, v1, nmlz);
      Vec cc0; fill(c, ci0, cc0, nmlz);
      Vec cc1; fill(c, ci1, cc1, nmlz);

      { // One of cc0 and cc1 is cc.
        Vec tmp;
        geo::axpbyz(1, cc0, -1, cc, tmp);
        Real mag = geo::norm(tmp);
        geo::axpbyz(1, cc1, -1, cc, tmp);
        mag = std::min(mag, geo::norm(tmp));
        if (mag > 1e-15) {
          // This is fatal, so report and immediately return.
          printf("neither cc0 nor cc1 is cc: %1.2e\n", mag);
          return;
        }
      }

      if (0) {
        // Edge center compared with the center computed from the
        // edge's vertices. This is not how the edge center is
        // computed.
        Vec ec1, tmp, vdif;
        geo::axpbyz(0.5, v0, 0.5, v1, ec1);
        geo::axpbyz(1, v0, -1, v1, vdif);
        geo::normalize(ec1);
        geo::axpbyz(1, ec, -1, ec1, tmp);
        ec_from_v_s(geo::norm(tmp)/geo::norm(vdif));
      }

      { // Edge center compared with the center computed from the
        // adjacent cell centers. This is the correct method.
        Vec ec1, tmp, vdif;
        geo::intersect_from_2linesegs(cc0, cc1, v0, v1, ec1);
        geo::axpbyz(1, ec, -1, ec1, tmp);
        geo::axpbyz(1, v0, -1, v1, vdif);
        ec_from_c_s(geo::norm(tmp)/geo::norm(vdif));
      }

      { // Check dvEdge (vtx-to-vtx distance). It's supposed to be
        // spherical length.
        Vec tmp;
        geo::axpbyz(1, v1, -1, v0, tmp);
        // Euclidean length
        const Real elen = R*geo::norm(tmp);
        const Real ere = std::abs(elen - g.dvEdge[ei])/g.dvEdge[ei];
        // Spherical length
        const Real slen = R*geo::calc_arc_length(v0, v1);
        const Real sre = std::abs(slen - g.dvEdge[ei])/g.dvEdge[ei];
        dvEdge_s(sre);
        if (sre > 1e-12) {
          if (out < 10) {
            printf("dvEdge: sre %1.2e ere %1.2e\n", sre, ere);
            out++;
          }
        }
      }

      { // Check dcEdge (cell-to-cell distance).
        const Real slen = R*geo::calc_arc_length(cc0, cc1);
        const Real sre = std::abs(slen - g.dcEdge[ei])/g.dcEdge[ei];
        dcEdge_s(sre);
      }

      { // Some angles.
        Vec cdif, vdif, tmp;
        geo::axpbyz(1, cc1, -1, cc0, cdif);
        geo::axpbyz(1, v1, -1, v0, vdif);
        geo::cross(cdif, vdif, tmp);
        { // Orthogonality of edge and two adjacent cell centers.
          const Real deg =
            std::abs(90 - (180.0/M_PI)*
                     std::asin(std::min(1.0,
                                        geo::norm(tmp)/
                                        (geo::norm(cdif)*geo::norm(vdif)))));
          ortho_s(deg);
        }
        { // Check angleEdge, which has range [-180, 180] M_PI/180.
          Vec axis = {0,0,1}, east, north, nml;
          Vec ec2; ec2[0] = ec[0]; ec2[1] = ec[1]; ec2[2] = 0;
          geo::cross(axis, ec2, east);
          geo::normalize(east);
          geo::cross(ec, east, north);
          geo::normalize(north);
          geo::cross(vdif, ec, nml);
          geo::normalize(nml);
          const Real theta = std::atan2(geo::dot(nml, north), geo::dot(nml, east));
          const Real deg = (180.0/M_PI)*std::min(
            std::min(std::abs(theta - g.angleEdge[ei]),
                     std::abs(theta - 2*M_PI - g.angleEdge[ei])),
            /**/     std::abs(theta + 2*M_PI - g.angleEdge[ei]));
          angle_s(deg);
        }
      }
    }

    // Vertex points on the sphere.
    for (int vi = 0; vi < ne; ++vi) {
      Vec v0; fill(v, t.get_c2v(ci,vi), v0);
      v_rad_s(std::abs(geo::norm(v0) - R)/R);
      geo::normalize(v0);
    }

    { // Check areaCell.
      Real vs[3*10];
      for (int vi = 0; vi < ne; ++vi)
        fill(v, t.get_c2v(ci,vi), vs + 3*vi, nmlz);
      const Real area = geo::calc_area_formula(vs, ne);
      area_s(std::abs(square(R)*area - g.areaCell[ci])/g.areaCell[ci]);
    }

    { // Check cell centroid.
      Vec cc1; calc_cell_centroid(m, ci, cc1);
      Vec diff; geo::axpbyz(1, cc, -1, cc1, diff);
      centroid_s(geo::norm(diff));
    }
  }

  printf(">> checks\n");
  c_rad_s.print("cc on sphere");
  e_rad_s.print("ec on sphere");
  v_rad_s.print("v on sphere");
  //ec_from_v_s.print("ec from v");
  ec_from_c_s.print("ec from c");
  dvEdge_s.print("dvEdge");
  dcEdge_s.print("dcEdge");
  area_s.print("area");
  angle_s.print("angle [deg]");
  ortho_s.print("ortho [deg]");
  centroid_s.print("centroid");
}

void get_angleEdge_error (const Topo& t, const CellCtrs& c, const EdgeCtrs& e,
                          const Vtxs& v, const MpasGeometry& g, Real* err) {
  Vec axis = {0,0,1};
  for (int ei = 0; ei < t.nedge; ++ei) {
    const int vi0 = t.get_e2v(ei,0), vi1 = t.get_e2v(ei,1);
    Vec ec; fill(e, ei, ec, true);
    Vec east; geo::cross(axis, ec, east); geo::normalize(east);
    Vec north; geo::cross(ec, east, north); geo::normalize(north);
    Vec v0; fill(v, vi0, v0, true);
    Vec v1; fill(v, vi1, v1, true);
    Vec vdif; geo::axpbyz(1, v1, -1, v0, vdif);
    Vec nml; geo::cross(vdif, ec, nml); geo::normalize(nml);
    const Real theta = std::atan2(geo::dot(nml, north), geo::dot(nml, east));
    const Real deg = (180.0/M_PI)*std::min(
      std::min(std::abs(theta - g.angleEdge[ei]),
               std::abs(theta - 2*M_PI - g.angleEdge[ei])),
      /**/     std::abs(theta + 2*M_PI - g.angleEdge[ei]));
    err[ei] = deg;
  }
}
