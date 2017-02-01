from cppstuff_to_expose import getLLH, getEvent
#import spectrum
from scipy.optimize import minimize
from scipy.optimize import brent
from scipy.optimize import minimize_scalar
from scipy.optimize import fmin
import numpy as np
from datetime import datetime

RESE = 200;
REST = 300;

triggEff = True
energyRes = True
mass = 1.0
distance = 5.0
events = 10
noise = pow(10,-3)*(10.0/REST)

# calculate llh for certain configuration
def llh(massi):
    return getLLH( float(massi), distance, events, triggEff, energyRes, noise, eventTime, eventEnergy)

# arrays to store the events
#eventEnergy = spectrum.intArray(int(events))
#eventTime = spectrum.intArray(int(events))
eventEnergy = np.zeros(int(events)).astype("int32")
eventTime = np.zeros(int(events)).astype("int32")

aa = datetime.now()
for i in range(1,10):
    # read in event
    getEvent(eventEnergy, eventTime, mass, distance, events, i)
    # minimize
    x_min = minimize_scalar(llh, bounds=(0,mass+1.0), method='bounded', options={'disp': True, 'xatol':0.005})
    print x_min.nfev, x_min.x, x_min.success, x_min.message
    #with open("DATA/masses_found_"+str(distance)+"Mpc_"+str(events)+"Events_"+str(mass)+"eV.txt", "a") as myfile:
    #    myfile.write(str(x_min.x) + '\n')

bb = datetime.now()
print "TIME", bb-aa

