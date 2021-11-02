# -*- coding: utf-8 -*-
"""
@author: alexa
"""


# Import data
from ddm import Sample
import pandas
import numpy as np
import pandas as pd

T_dur = 10 # max duration of trial
Ntrial = 60 # number of trials per condition
Nsimul = 500 # number of simulations
DriftA = 0.8 #drift in condition A
DriftB= 0.7 #drift in condition B
nondectime = .3 #non-decision time (s)

from ddm import Model, Fittable, Drift
from ddm.functions import fit_adjust_model
from ddm.models import NoiseConstant, BoundConstant, OverlayNonDecision

class DriftCondition(Drift):
    name = "Drift depends on condition"
    required_parameters = ["driftA", "driftB"] # <-- Drift parameters for each condition
    required_conditions = ["condition"] # <-- Task parameters ("conditions"). Should be the same name as in the sample.
    
    # We must always define the get_drift function, which is used to compute the instantaneous value of drift.
    def get_drift(self, conditions, **kwargs):
        return self.driftA if conditions['condition']==1 else self.driftB


model = Model(name='Simple model',
             drift=DriftCondition(driftA=DriftA, driftB = DriftB),
              noise=NoiseConstant(noise=1),
              bound=BoundConstant(B=1),
              overlay=OverlayNonDecision(nondectime=nondectime),
              dx=.001, dt=.01, T_dur=T_dur)

model_fit = Model(name='Simple model (fitted)',
                  drift=DriftCondition(driftA=Fittable(minval=.3, maxval=3), driftB=Fittable(minval=.3, maxval=3)),
                  noise=NoiseConstant(noise=1),
                  bound=BoundConstant(B=Fittable(minval=.5, maxval = 3)),
                  overlay=OverlayNonDecision(nondectime=Fittable(minval=0, maxval=1)),
                  dx=.001, dt=.01, T_dur=T_dur)

drift_diff = np.empty(Nsimul,)

# simulate and fit
for i in range(Nsimul):
    
    #simulate trials
    solA = model.solve(conditions ={"condition" : 1})
    sampA = solA.resample(Ntrial)
    dfA = sampA.to_pandas_dataframe()
    solB = model.solve(conditions ={"condition" : 0})
    sampB = solB.resample(Ntrial)
    dfB = sampB.to_pandas_dataframe()
    
    df_all = pd.concat([dfA, dfB], axis=0) # concatenate both
    sample = Sample.from_pandas_dataframe(df_all, rt_column_name="RT", correct_column_name="correct")    
        
    # Fit
    fit_model = fit_adjust_model(sample=sample, model=model_fit, verbose=False)
    par_hat = fit_model.get_model_parameters()
    drift_diff[i] = float(par_hat[0]) - float(par_hat[1]) # estimated difference in drift

    if (i%10==0):
        print('o')

# export as CSV
pd.DataFrame(drift_diff).to_csv("DDMeffectsize_simulations.csv")