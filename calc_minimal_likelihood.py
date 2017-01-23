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
filenumber = 1
noise = pow(10,-5);

# calculate llh for certain configuration
def llh(massi):
    return spectrum.getLLH( float(massi), distance, events, triggerEfficiency, energyRes, noise, eventTime, eventEnergy)

# test function
def x2(x):
    return (x-2)*(x-2)+3;

# get trigger efficiency
triggerEfficiency = spectrum.doubleArray(601)
spectrum.getTriggerEfficiency(triggerEfficiency, triggEff)

# read in event - TODO: loop over all events
eventEnergy = spectrum.intArray(160)
eventTime = spectrum.intArray(160)
spectrum.getEvent(eventEnergy, eventTime, mass, distance, events, 2);

'''
a = datetime.now()
test2 = minimize(llh, 1.0, method='nelder-mead', options={'xtol': 0.01, 'disp': True})
print 'result 2', test2.x[0]
b = datetime.now()
print "TIME 2:", b-a

a = datetime.now()
xlow,path = fmin(llh,1.0,retall=True, xtol=0.01,)
print xlow,path
b = datetime.now()

print "TIME 3:", b-a
'''

# seems to be the best so far, need to set bounds
aa = datetime.now()
x_min = minimize_scalar(llh, bounds=(0,mass+1.0), method='bounded', options={'disp':1,'xatol':0.005})
print x_min
bb = datetime.now()
print "TIME", bb-aa
