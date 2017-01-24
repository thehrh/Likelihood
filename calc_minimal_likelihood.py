import spectrum
from scipy.optimize import minimize
from scipy.optimize import brent
from scipy.optimize import minimize_scalar
from scipy.optimize import fmin
import numpy as np
from datetime import datetime

triggEff = True
energyRes = True
mass = 1.0
distance = 1.0
events = 160
noise = pow(10,-5);

# calculate llh for certain configuration
def llh(massi):
    return spectrum.getLLH( float(massi), distance, events, triggEff, energyRes, noise, eventTime, eventEnergy)

# arrays to store the events
eventEnergy = spectrum.intArray(160)
eventTime = spectrum.intArray(160)

aa = datetime.now()
for i in range(1,100):
    # read in event
    spectrum.getEvent(eventEnergy, eventTime, mass, distance, events, i);
    # minimize
    x_min = minimize_scalar(llh, bounds=(0,mass+1.0), method='bounded', options={'disp':1,'xatol':0.005})
    print x_min.x
    with open("TEST_minimizer.txt", "a") as myfile:
        myfile.write(str(x_min.x) + '\n')

bb = datetime.now()
print "TIME", bb-aa

