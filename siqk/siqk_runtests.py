#!/usr/bin/python

import os

quick = True

stride = 1
biggest = 1111

xlates = [4.2*10**f for f in range(-17, 0, stride)]
xlates.append(0)

ylates = [0]

angles = xlates

fails = []
cnt = 0

# Test 1
for n in [4, 20, 40, 79]:
    if quick and n > 20: break
    for angle in angles:
        cmd = ('OMP_NUM_THREADS=8 ./siqk_test --testno 1 --angle {angle:1.15e} -n {n:d}'.
               format(angle=angle, n=n))
        stat = os.system(cmd + ' |& grep PASSED &> /dev/null')
        if stat:
            fails.append(cmd)
        else:
            cnt += 1
    print len(fails)

# Test 0
for n in [4, 50, 511, biggest]:
    if quick and n > 50: break
    for angle in angles:
        for xlate in xlates:
            for ylate in ylates:
                cmd = ('OMP_NUM_THREADS=8 ./siqk_test --testno 0 --xlate {xlate:1.15e} --ylate {ylate:1.14e} --angle {angle:1.15e} -n {n:d}'.
                       format(xlate=xlate, ylate=ylate, angle=angle, n=n))
                stat = os.system(cmd + ' |& grep PASSED &> /dev/null')
                if stat:
                    fails.append(cmd)
                else:
                    cnt += 1
    print len(fails)
    
if len(fails) > 0:
    print 'FAILED'
    for f in fails:
        print f
else:
    print 'PASSED ({0:d})'.format(cnt)
