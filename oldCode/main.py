import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

global N, A
A=1
Narr = np.arange(20000, 30001, 10000)
print(len(Narr))
for iN in range(0, len(Narr)):
    N = Narr[iN]
    # execfile('settings.py')