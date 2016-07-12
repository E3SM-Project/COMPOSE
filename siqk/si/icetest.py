#!/usr/bin/python

import os

stride = 1

xlates = [4.2*10**f for f in range(-17, 0, stride)]
xlates.append(0)

ylates = [0]

angles = xlates

fails = []
cnt = 0

for n in [21]:
    for angle in angles:
        for xlate in xlates:
            for ylate in ylates:
                cmd = ('./test.exe --plane --xlate {xlate:1.15e} --ylate {ylate:1.14e} --angle {angle:1.15e} -n {n:d}'.
                       format(xlate=xlate, ylate=ylate, angle=angle, n=n))
                stat = os.system(cmd + ' |& grep PASSED &> /dev/null')
                if stat:
                    fails.append(cmd)
                else:
                    cnt += 1
    print len(fails)
