import mobcal
import numpy as np

RngGen=np.random.RandomState(seed=43)
mobcal.xrand=RngGen.rand
mobcal.mobcal()