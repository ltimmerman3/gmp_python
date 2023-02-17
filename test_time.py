import time
from math import comb
from intderiv_assist import nCk

t0 = time.time()
comb(10,3)
t1 = time.time()
nCk(10,3)
t2 = time.time()
print('comb: {}\tnCk: {}'.format(t1-t0,t2-t1))
