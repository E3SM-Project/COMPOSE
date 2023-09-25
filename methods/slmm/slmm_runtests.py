#!/usr/bin/python

import os, sys, re, optparse

def readall (fn):
    # Shorthand for reading in all the text in a file.
    try:
        with open(fn, 'r') as f:
            text = f.read()
    except:
        text = ''
    return text

def writeall (text, fn, for_real):
    if for_real:
        with open(fn, 'w') as f:
            f.write(text)

def parse_one_liner (text):
    class struct:
        pass
    hits = re.findall('<OL>.*', text)
    hits = re.findall('l2 (?P<l2>[^ ]*) .* cv re (?P<cv>[^ ]*)' +
                      '.* cvgll re (?P<cvgll>[^ ]*)' +
                      '.* mo min [0-9.e\-+]+ (?P<min>[^ ]*) .* ' +
                      'max [0-9.e\-+]+ (?P<max>[^ ]*)', hits[0])
    o = struct
    o.l2 = float(hits[0][0])
    o.cv = float(hits[0][1])
    o.cv_gll = float(hits[0][2])
    o.mo_min = float(hits[0][3])
    o.mo_max = float(hits[0][4])
    return o

def runtest (cmd):
    outfn = 'runtests.tmp'
    os.system(cmd + ' > ' + outfn + ' 2>&1')
    return readall(outfn)

long_output = False

def print_test (cmd):
    print_test.ctr += 1
    ll = 87;
    if long_output:
        ll = 240
    if len(cmd) > ll:
        cmd = cmd[len(cmd)-ll+1:]
    fmt = '{{0:.<{0:d}s}}'.format(ll+1)
    print('{:3d} '.format(print_test.ctr) + fmt.format(cmd + ' '), end='')
    sys.stdout.flush()
print_test.ctr = 0;

def print_result (passed):
    if not passed:
        print('***FAILED')
        return 1
    else:
        print('   PASSED')
        return 0

def check_passed (cmd):
    print_test(cmd)
    out = runtest(cmd)
    hits = re.findall('PASSED', out)
    passed = len(hits) > 0
    return print_result(passed)

def check_errs (cmd, l2_err, cv=10, cv_gll=10, min=-float('Inf'), max=float('Inf'),
                l2_err_is_0=False):
    print_test(cmd)
    out = runtest(cmd)
    o = parse_one_liner(out)
    passed = ((o.l2 > 0 or l2_err_is_0) and o.l2 <= l2_err
              and o.cv <= cv
              and o.cv_gll <= cv_gll
              and o.mo_min >= min and o.mo_max <= max)
    result = print_result(passed)
    if not passed:
        print('   ' + cmd)
        print(('    l2 {:1.2e} cv {:1.2e} cv_gll {:1.2e} mo_min {:1.2e} mo_max {:1.2e}' +
               ' but l2_err {:1.2e} cv {:1.2e} cv_gll {:1.2e} min {:1.2e} max {:1.2e}').
            format(o.l2, o.cv, o.cv_gll, o.mo_min, o.mo_max,
                   l2_err, cv, cv_gll, min, max))
    return result

p = optparse.OptionParser()
p.add_option('-l', '--long', dest='long', action='store_true', default=False,
             help='Long-line output.')
opts, args = p.parse_args()
long_output = opts.long

try: os.mkdir('tmp')
except: pass

nerr = 0
# Unit tests.
nerr += check_passed('./slmm_test -q -c test_make_cubedsphere')
nerr += check_passed('./slmm_test -q -c test_gll')
nerr += check_passed('./slmm_test -q -c test_gll_2d')
nerr += check_passed('./slmm_test -q -c test_time_int')
nerr += check_passed('./slmm_test -q -c test_make_gll_mesh')
nerr += check_passed('./slmm_test -q -c test_make_gll_subcell_mesh')
nerr += check_passed('./slmm_test -q -c test_qp_limiter')
nerr += check_passed('./slmm_test -q -c test_face_tree')
nerr += check_passed('./slmm_test -q -c test_spf')
nerr += check_passed('./slmm_test -q -c test_nla')
nerr += check_passed('./slmm_test -q -c test_mass_matrix')
#nerr += check_passed('./slmm_test -q -c test_fit_extremum')

# Test classical semi-Lagrangian with global filters QLT, CAAS, min-norm2.
base = ('./slmmir -method {method:s} -ode divergent -ic slottedcylinders ' +
        '-ic cosinebells -ic gaussianhills -we 0 -np {np:d} -dmc f -mono {mono:s} ' +
        '-nsteps 12 -ne {ne:d}')
nerr += check_errs(base.format(method='pisl', np=4, ne=10, mono='qlt'),
                   3.34e-1, cv_gll=5e-14, min=0.1, max=1)  # rho is also done with ISL
nerr += check_errs(base.format(method='pisl', np=6, ne=6, mono='qlt'),
                   3.34e-1, cv_gll=5e-14, min=0.1, max=1)  # rho is also done with ISL
nerr += check_errs(base.format(method='isl', np=4, ne=10, mono='qlt'),
                   3.47e-1, cv_gll=5e-14, min=0.1, max=1)  # rho is remapped
nerr += check_errs(base.format(method='pisl', np=4, ne=10, mono='qlt-pve'),
                   3.36e-1, cv_gll=5e-14, min=0, max=2)    # >= 0 constraint only
nerr += check_errs(base.format(method='pisl', np=4, ne=10, mono='caas'),
                   3.47e-1, cv_gll=5e-14, min=0.1, max=1)  # rho is also done with ISL
nerr += check_errs(base.format(method='isl', np=4, ne=10, mono='caas'),
                   3.47e-1, cv_gll=5e-14, min=0.1, max=1)  # rho is remapped
nerr += check_errs(base.format(method='isl', np=4, ne=10, mono='mn2'),
                   3.47e-1, cv_gll=5e-14, min=0.1, max=1)  # rho is remapped
#   Tracer consistency test. Apply ISL to constant q but remap rho.
nerr += check_errs('./slmmir -method isl -ode divergent -ic constant -we 0 -np 4 ' +
                   '-dmc f -mono qlt -rit -nsteps 12 -ne 10',
                   3e-15, cv_gll=1e-13, min=0.42, max=0.42, l2_err_is_0=True)

# Test ISL with p-refinement.
base = ('./slmmir -method pisl -ode divergent -ic gaussianhills ' +
        '-we 0 -np {np:d} -dmc f -mono {mono:s} ' +
        '-nsteps 12 -ne {ne:d} -timeint {timeint:s}')
nerr += check_errs(base.format(np=12, ne=3, mono='none', timeint='interp'),
                   9.939e-3)
nerr += check_errs(base.format(np=12, ne=3, mono='none', timeint='exact'),
                   8.793e-3)
base = ('./slmmir -method pisl -ode divergent -ic slottedcylinders '
        '-we 0 -np {np:d} -dmc f -mono {mono:s} ' +
        '-nsteps 12 -ne {ne:d} -timeint interp')
nerr += check_errs(base.format(np=12, ne=3, mono='caas', timeint='interp'),
                   2.896e-1, cv_gll=5e-14, min=0.1, max=1)

# ISL with p-refinement and separate t and v meshes.
base = ('./slmmir -method pisl -ode divergent -ic gaussianhills ' +
        '-we 0 -rit -dmc {dmc:s} -mono {mono:s} -lim {lim:s} -nsteps 13 -T 12 ' +
        '-ne 6 -np 8 -timeint interp -prefine {prefine:d} -d2c')
nerr += check_errs(base.format(prefine=0, dmc='es', mono='caas', lim='caas'),  5.968e-03, cv=2e-14)
nerr += check_errs(base.format(prefine=5, dmc='es', mono='caas', lim='caas'),  5.885e-03, cv=4e-14)
nerr += check_errs(base.format(prefine=0, dmc='eh', mono='caas', lim='caas'),  5.968e-03, cv_gll=2e-14)
nerr += check_errs(base.format(prefine=5, dmc='eh', mono='caas', lim='caas'),  5.886e-03, cv_gll=2e-14)
# new global-only method
nerr += check_errs(base.format(prefine=0, dmc='es', mono='caas-node', lim='caas'),  5.968e-03, cv=2e-14)
nerr += check_errs(base.format(prefine=5, dmc='es', mono='caas-node', lim='caas'),  5.885e-03, cv=4e-14)
nerr += check_errs(base.format(prefine=0, dmc='eh', mono='caas-node', lim='caas'),  5.968e-03, cv_gll=2e-14)
nerr += check_errs(base.format(prefine=5, dmc='eh', mono='caas-node', lim='caas'),  5.886e-03, cv_gll=2e-14)
# don't break the no prop preserve case
nerr += check_errs(base.format(prefine=5, dmc='es', mono='none', lim='none'),  4.2e-03)
# rotate the grid
rbase = base + ' -rotate-grid'
nerr += check_errs(rbase.format(prefine=5, dmc='eh', mono='caas-node', lim='caas'),  5.886e-03, cv_gll=2e-14)
# GllOffsetNodal
base += ' -basis GllOffsetNodal'
nerr += check_errs(base.format(prefine=5, dmc='es', mono='caas', lim='caas'),  5.885e-03, cv=4e-14)
nerr += check_errs(base.format(prefine=5, dmc='eh', mono='caas', lim='caas'),  5.886e-03, cv_gll=2e-14)
nerr += check_errs(base.format(prefine=5, dmc='es', mono='caas-node', lim='caas'),  5.885e-03, cv=4e-14)
nerr += check_errs(base.format(prefine=5, dmc='eh', mono='caas-node', lim='caas'),  5.886e-03, cv_gll=2e-14)

base = './slmmir -nsteps 12 -ne 10 -we 0 -ode divergent -ic gaussianhills '

# DSS for QOF rho, ISL tracer, with QLT.
nerr += check_errs(base + '-np 3 -d2c -method isl -dmc f -mono qlt',
                   9.05e-2, cv_gll=2e-14)

# Cell-integrated method basics.
nerr += check_errs(base + '-np 3', 2.43e-2, 1e-14)
nerr += check_errs(base + '-np 3 -xyz -mono qlt',
                   3.18e-2, 4e-15, min=1.495e-08, max=9.518e-01)
nerr += check_errs(base + '-np 3 -xyz -mono caas',
                   3.18e-2, 4e-15, min=1.495e-08, max=9.518e-01)
nerr += check_errs(base + '-np 3 -xyz -mono mn2',
                   3.18e-2, 4e-15, min=1.495e-08, max=9.518e-01)
nerr += check_errs(base + '-np 3 -xyz -d2c', 3.64e-2, 3e-15)
nerr += check_errs(base + '-np 4 -xyz -d2c', 1.02e-2, 8e-15)
nerr += check_errs(base + '-np 4 -xyz -d2c -method cdg', 1.02e-2, 3.5e-15)

# Limiter.
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -we 0 -ode divergent ' +
                   '-ic slottedcylinders -np 4 -mono qlt -method ir',
                   3.0e-1, cv=3e-14, min=0.1, max=1.0)
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -we 0 -ode divergent ' +
                   '-ic slottedcylinders -np 4 -mono qlt -method ir -lim caas',
                   3.0e-1, cv=3e-14, min=0.1, max=1.0)
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -we 0 -ode divergent ' +
                   '-ic slottedcylinders -np 4 -mono qlt -method cdg',
                   3.03e-1, cv=3e-14, min=0.1, max=1.0)
# Multiple tracers.
nerr += check_errs(base + '-np 4 -ic correlatedcosinebells 2', 1.02e-2, 2e-7)
# Local DMC with internal mass definition.
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -we 0 -ode divergent ' +
                   '-ic gaussianhills -np 4 -dmc es -method ir',
                   9.1e-3, cv=2e-13)
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -we 0 -ode divergent ' +
                   '-ic gaussianhills -np 4 -dmc es -method cdg',
                   9.1e-3, cv=2e-13)
# Local DMC with Homme mass definition.
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -we 0 -ode divergent ' +
                   '-ic gaussianhills -np 4 -dmc eh',
                   9.1e-3, cv_gll=5e-15)
# Global (weaker than local) DMC with Homme mass definition.
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -we 0 -ode divergent ' +
                   '-ic gaussianhills -np 4 -dmc geh',
                   9.1e-3, cv_gll=2e-14)
# Local DMC, limiter, internal mass def.
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -we 0 -ode divergent ' +
                   '-ic slottedcylinders -np 4 -mono qlt -dmc es',
                   3.1e-1, cv=2.3e-13, min=0.1, max=1.0)
# Local DMC, limiter, Homme mass def.
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -we 0 -ode divergent ' +
                   '-ic slottedcylinders -np 4 -mono qlt -dmc eh',
                   3e-1, cv_gll=5e-14, min=0.1, max=1.0)
# Local DMC, facet transport.
nerr += check_errs(base + '-np 4 -dmc f', 1.42e-2, cv_gll=6e-14)
nerr += check_errs('./slmmir -nsteps 12 -ne 30 -we 0 -ode divergent ' +
                   '-ic gaussianhills -np 2 -dmc f',
                   6.49e-2, cv_gll=1.4e-13)
# With limiter.
nerr += check_errs('./slmmir -nsteps 96 -ne 5 -we 0 -ode divergent -method cdg ' +
                   '-ic slottedcylinders -np 4 -mono qlt -dmc f',
                   4.6e-1, cv_gll=4e-14, min=0.1, max=1.0)
nerr += check_errs('./slmmir -nsteps 96 -ne 5 -we 0 -ode divergent -method cdg ' +
                   '-ic slottedcylinders -np 4 -mono qlt -dmc f -lim caas',
                   4.6e-1, cv_gll=4e-14, min=0.1, max=1.0)
nerr += check_errs('./slmmir -nsteps 96 -ne 5 -we 0 -ode divergent -method cdg ' +
                   '-ic slottedcylinders -np 4 -mono qlt -dmc f -lim caags',
                   4.6e-1, cv_gll=4e-14, min=0.1, max=1.0)
nerr += check_errs('./slmmir -nsteps 96 -ne 5 -we 0 -ode divergent -method ir ' +
                   '-ic slottedcylinders -np 4 -mono qlt -dmc f',
                   4.6e-1, cv_gll=4e-14, min=0.1, max=1.0)
# Add an equality constraint to nail DMC even more. In addition, output scalar
# measurements by time step (this is just a test that it runs).
nerr += check_errs('./slmmir -nsteps 96 -ne 5 -we 0 -ode divergent -method cdg ' +
                   '-ic slottedcylinders -np 4 -mono qlt -dmc ef -o rittest -rit',
                   4.6e-1, cv_gll=2e-14, min=0.1, max=1.0)
nerr += check_errs('./slmmir -nsteps 96 -ne 5 -we 0 -ode divergent -method ir ' +
                   '-ic slottedcylinders -np 4 -mono qlt -dmc ef -o rittest -rit',
                   4.6e-1, cv_gll=2e-14, min=0.1, max=1.0)
nerr += check_errs('./slmmir -nsteps 96 -ne 15 -we 0 -ode divergent ' +
                   '-ic slottedcylinders -np 2 -mono qlt -dmc ef -o rittest -rit',
                   4.5e-1, cv_gll=2.2e-14, min=0.1, max=1.0)
# Test the more complicated mono method.
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -we 0 -ode divergent ' +
                   '-ic gaussianhills -ic slottedcylinders -np 4 -mono qlt -dmc f',
                   1.5e-2, cv_gll=8e-14, min=0, max=0.957)
# 3-1 subcell mesh, with new vertices at GLL points.
nerr += check_errs('./slmmir -nsteps 96 -ne 5 -we 0 -ode divergent -tq 4 ' +
                   '-ic slottedcylinders -np 4 -mesh gllsubcell -mono qlt -dmc ef',
                   4.6e-1, cv_gll=2e-14, min=0.1, max=1.0)
# 3-1 subcell mesh, with new vertices at non-GLL points.
nerr += check_errs('./slmmir -nsteps 96 -ne 5 -we 0 -ode divergent -tq 4 ' +
                   '-ic slottedcylinders -np 4 -mesh runisubcell -mono qlt -dmc ef',
                   4.5e-1, cv_gll=2e-14, min=0.1, max=1.0)
# Same, but now looking for accuracy difference.
nerr += check_errs('./slmmir -nsteps 12 -ne 5 -we 0 -ode divergent -tq 4 ' +
                   '-ic gaussianhills -np 4 -mesh gllsubcell -mono qlt -dmc ef',
                   7.40e-2, cv_gll=9e-15, min=0, max=0.96)
nerr += check_errs('./slmmir -nsteps 12 -ne 5 -we 0 -ode divergent -tq 4 ' +
                   '-ic gaussianhills -np 4 -mesh runisubcell -mono qlt -dmc ef',
                   5.41e-2, cv_gll=5e-15, min=0, max=0.96)
# We can subdivide cells arbitrarily with runisubcell.
nerr += check_errs('./slmmir -nsteps 12 -ne 2 -we 0 -ode divergent -tq 4 ' +
                   '-ic gaussianhills -np 10 -mesh runisubcell -mono qlt -dmc ef',
                   3.5e-2, cv_gll=3e-15, min=0, max=0.96)
# Tracer-decoupled CMBC tests.
base = ('./slmmir -nsteps 12 -ne 10 -np 4 -ode divergent ' +
        '-ic gaussianhills -ic slottedcylinders -ic cosinebells ' +
        '-ic correlatedcosinebells -ic xyztrig -dmc {0:s} -mono {mono:s} -we 0')
#  This method also is intended to handle tracer consistency, but I haven't put
#  together a test for that yet. So test just CMBC.
nerr += check_errs(base.format('f', mono='qlt'), 1.45e-2, cv_gll=6e-14, min=1.495e-8, max=0.956)
nerr += check_errs(base.format('es', mono='qlt'), 9.18e-3, cv=2e-13, min=1.495e-8, max=0.956)
nerr += check_errs(base.format('eh', mono='qlt'), 9.18e-3, cv_gll=1e-14, min=1.495e-8, max=0.956)
#  Test that if rho is perturbed, a constant q stays a constant.
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -np 4 -ode nondivergent ' +
                   '-ic constant -dmc ef -mono qlt -we 0 --perturb-rho 0.05',
                   1e-6, cv_gll=5e-14, min=0.42, max=0.42, l2_err_is_0=True)
nerr += check_errs('./slmmir -nsteps 12 -ne 10 -np 4 -ode divergent ' +
                   '-ic constant -dmc ef -mono qlt -we 0 --perturb-rho 0.05',
                   1e-6, cv_gll=5e-14, min=0.42, max=0.42, l2_err_is_0=True)

print('{0:d} tests failed'.format(nerr))
