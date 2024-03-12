#include "slmmir_p_refine.hpp"
#include "slmmir_lauritzen_diag.hpp"

#include <vector>

extern "C" {
  void correlation_diag(const Real* f1, const Real* f2, int K, const Real* dA,
                        Real* real_mixing, Real* overshooting, Real* range_pres_unmixing);
  void filament_diag(int K, const Real* f1, const Real* dA, Real* fila_t0,
                     int linit, Real* thresholds, Real* fila_tf);
}

static void write_mixing_plot_data (const Real* cb, const Real* ccb, const Int n,
                                    const D2Cer& d2cer, const std::string& out_prefix) {
  const Int cnn = d2cer.get_cnn();
  std::vector<Real> cgdata(2*cnn);
  auto* const cbd = cgdata.data();
  auto* const ccbd = cgdata.data() + cnn;
  d2cer.d2c(cb, cbd);
  d2cer.d2c(ccb, ccbd);
  const auto fname = out_prefix + "-lauritzen-diag.dat";
  FILE* fid = fopen(fname.c_str(), "wb");
  if ( ! fid) {
    printf("Could not open %s; exiting.\n", fname.c_str());
    exit(-1);
  }
  std::vector<float> cgsp(2*cnn);
  for (int i = 0; i < 2*cnn; ++i) cgsp[i] = cgdata[i];
  for (int i = 0; i < 2*cnn; ++i) if (std::isnan(cgsp[i])) prc(i);
  auto* const cbdsp = cgsp.data();
  auto* const ccbdsp = cgsp.data() + cnn;
  fwrite(&cnn, sizeof(Int), 1, fid);
  fwrite(cbdsp, sizeof(*cbdsp), cnn, fid);
  fwrite(ccbdsp, sizeof(*ccbdsp), cnn, fid);
  fclose(fid);
}

struct LauritzenDiag::Impl {
  Int nsteps_per_12days, len;
  std::vector<Real> thresholds, fila_t0, fila_tf, fila_t0_me, fila_tf_me;
  Int tidx_cosbells, tidx_corrcosbells;
  bool expensive_io, ran_cor, ran_fil, ran_me;
  Real real_mixing, overshooting, range_pres_unmixing;
  Real real_mixing_me, overshooting_me, range_pres_unmixing_me;
};

LauritzenDiag
::LauritzenDiag (const Int nsteps_per_12days, const Int len,
                 const std::vector<gallery::InitialCondition::Shape>& ics,
                 const Real* tracer_data, const Real* dA,
                 const bool expensive_io) {
#ifdef SLMMIR_LAURITZEN_DIAG
  p = std::make_shared<Impl>();
  p->expensive_io = expensive_io;
  p->nsteps_per_12days = nsteps_per_12days;
  p->len = len;
  p->thresholds.resize(100);
  p->fila_t0.resize(100); p->fila_tf.resize(100);
  p->fila_t0_me.resize(100); p->fila_tf_me.resize(100);
  p->tidx_cosbells = p->tidx_corrcosbells = -1;
  p->ran_cor = p->ran_fil = p->ran_me = false;
  for (size_t i = 0; i < ics.size(); ++i) {
    if (ics[i] == gallery::InitialCondition::Shape::CosineBells)
      p->tidx_cosbells = i;
    else if (ics[i] == gallery::InitialCondition::Shape::CorrelatedCosineBells)
      p->tidx_corrcosbells = i;
  }
  if (p->tidx_cosbells >= 0)
    filament_diag(len, tracer_data + p->tidx_cosbells*len, dA,
                  p->fila_t0.data(), 1, nullptr, nullptr);
#endif
}

bool LauritzenDiag
::run (const Int step, const Real* tracer_data, const Real* dA,
       const D2Cer& d2cer, const std::string& out_prefix) {
#ifdef SLMMIR_LAURITZEN_DIAG
  if (p->tidx_cosbells < 0) return false;
  const bool day6 = ((p->nsteps_per_12days % 2 == 0) &&
                     (step + 1 == p->nsteps_per_12days/2));
  if ( ! day6) return false;
  filament_diag(p->len, tracer_data + p->tidx_cosbells * p->len, dA,
                p->fila_t0.data(), 0, p->thresholds.data(), p->fila_tf.data());
  p->ran_fil = true;
  if (p->tidx_corrcosbells < 0) return false;
  correlation_diag(tracer_data + p->tidx_cosbells * p->len,
                   tracer_data + p->tidx_corrcosbells * p->len, p->len, dA,
                   &p->real_mixing, &p->overshooting, &p->range_pres_unmixing);
  if (p->expensive_io)
    write_mixing_plot_data(tracer_data + p->tidx_cosbells * p->len,
                           tracer_data + p->tidx_corrcosbells * p->len,
                           p->len, d2cer, out_prefix);
  p->ran_cor = true;
  return p->ran_fil && p->ran_cor;
#endif
  return false;
}

void LauritzenDiag
::run_me (const Int n, const Real* tracer_data, const Real* dA, const bool first) {
#ifdef SLMMIR_LAURITZEN_DIAG
  if ( ! (p->ran_fil || first) || p->ran_me) return;
  if ( ! first) p->ran_me = true;
  if (first)
    filament_diag(n, tracer_data + p->tidx_cosbells * n, dA,
                  p->fila_t0_me.data(), 1, nullptr, nullptr);
  else
    filament_diag(n, tracer_data + p->tidx_cosbells * n, dA,
                  p->fila_t0_me.data(), 0, p->thresholds.data(),
                  p->fila_tf_me.data());
  if (p->tidx_corrcosbells < 0) return;
  correlation_diag(tracer_data + p->tidx_cosbells * n,
                   tracer_data + p->tidx_corrcosbells * n, n, dA,
                   &p->real_mixing_me, &p->overshooting_me,
                   &p->range_pres_unmixing_me);
#endif
}

void LauritzenDiag::print () {
  if (p->ran_cor) {
    printf("L    l_r %8.2e l_u %8.2e l_o %8.2e\n",
           p->real_mixing, p->range_pres_unmixing, p->overshooting);    
    if (p->ran_me)
      printf("L me l_r %8.2e l_u %8.2e l_o %8.2e\n",
             p->real_mixing_me, p->range_pres_unmixing_me, p->overshooting_me);
  }
  if (p->ran_fil) {
    printf("L    thr");
    Int n = 0;
    for (size_t i = 0; i < p->thresholds.size(); ++i) {
      if (p->thresholds[i] < 0) {
        n = i;
        break;
      }
      printf(" %1.3f", p->thresholds[i]);
    }
    printf("\nL    fil");
    for (Int i = 0; i < n; ++i) printf(" %1.2f", p->fila_tf[i]);
    printf("\n");
    if (p->ran_me) {
      printf("\nL me fil");
      for (Int i = 0; i < n; ++i) printf(" %1.2f", p->fila_tf_me[i]);
      printf("\n");      
    }
  }
}
