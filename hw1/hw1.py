# Activate any packages and imports
from acoli import acoli_hist
import numpy as np
# Package for 3-D
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Heat map for food function 1
rttp = 0.5
ttrp = 0.5
trials = 1
steps = 50
bins = 100
acoli_hist([0,0], 2, rttp, ttrp, trials, steps, bins, True)



