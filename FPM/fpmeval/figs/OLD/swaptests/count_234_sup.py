#!/usr/bin/env python

import sys, os
CNT = {}
for f in sys.argv[1:]:
    cmd = "grep 2 %s | grep 3 | grep 4 | wc" % f
    res = os.popen(cmd).readlines()[0].split()[0]
    if res not in CNT: CNT[res] =0
    CNT[res] += 1

print CNT
