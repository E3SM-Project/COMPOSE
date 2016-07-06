#!/usr/bin/python

import os, sys, re

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
    hits = re.findall('l2 (?P<l2>[^ ]*) .* cv re (?P<cv>[^ ]*) ', hits[0])
    o = struct
    o.l2 = float(hits[0][0])
    o.cv = float(hits[0][1])
    return o

def runtest (cmd):
    outfn = 'runtests.tmp'
    os.system(cmd + ' &> ' + outfn)
    return readall(outfn)

def print_test (cmd):
    print '{0:.<70s}'.format(cmd + ' '),

def print_result (passed):
    if not passed:
        print '***FAILED'
        return 1
    else:
        print '   PASSED'
        return 0

def check_passed (cmd):
    print_test(cmd)
    out = runtest(cmd)
    hits = re.findall('PASSED', out)
    passed = len(hits) > 0
    return print_result(passed)

def check_errs (cmd, l2_err, cv=10):
    print_test(cmd)
    out = runtest(cmd)
    o = parse_one_liner(out)
    passed = o.l2 <= l2_err and o.cv <= cv
    return print_result(passed)

nerr = 0
nerr += check_passed('./slmm_test -q -c test_make_cubedsphere')
nerr += check_passed('./slmm_test -q -c test_gll')
nerr += check_passed('./slmm_test -q -c test_time_int')
nerr += check_passed('./slmm_test -q -c test_make_gll_mesh')

base = './slmmir -nsteps 12 -ne 10 -we 0 -ode divergent -ic gaussianhills '
nerr += check_errs(base + '-np 3', 6.3e-3, 1e-5)
nerr += check_errs(base + '-np 3 -xyz', 6.3e-3, 1e-5)
nerr += check_errs(base + '-np 3 -xyz -d2c', 8.8e-3, 1e-5)
nerr += check_errs(base + '-np 4 -xyz -d2c', 5e-3, 2e-7)

print '{0:d} tests failed'.format(nerr)
