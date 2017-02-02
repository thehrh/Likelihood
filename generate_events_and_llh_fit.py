from cppstuff_to_expose import getLLH, getEvent, createSpectrum
#import spectrum
from scipy.optimize import minimize
from scipy.optimize import brent
from scipy.optimize import minimize_scalar
from scipy.optimize import fmin
import numpy as np
from datetime import datetime
import argparse


# make binning configurable?
RESE = 200;
REST = 300;

#triggEff = True
#energyRes = True
#mass = 1.0
#distance = 5.0
#events = 10
noise = pow(10,-3)*(10.0/REST)

# calculate llh for certain configuration
def llh(massi):
    return getLLH( float(massi), distance, nevts, useTriggerEff, useEnergyRes, noise, eventTime, eventEnergy)

def sample_spectrum():
    nevts_generated = 0
    evts = []
    while nevts_generated < nevts:
	# seed: make reproducible and non-redundant
        randE = np.random.randint(0, RESE-1)
        randT = np.random.randint(0, REST)
        rand_check = np.random.rand()*spectrum_max
        if spectrum[randT*(RESE-1)+randE] >= rand_check:
            evts.append([randE, randT])
            nevts_generated += 1
    return evts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
		description="Simulate SN Spectrum and fit pseudo-experiments.")
    parser.add_argument("-m", "--mass", default=1.0,
			help="Neutrino mass to simulate")
    parser.add_argument("-d", "--distance", default=5.0,
			help="SN Distance (in Mpc) to simulate")
    parser.add_argument("-n", "--nevents", required=True, type=float,
			help="No. of events. IMPROVE ME!")
    parser.add_argument("--perfect-trigger", dest='triggEff',
			default=True, action='store_false',
			help="Assume fully  eff. trigger across all energies.")
    parser.add_argument("--perfect-reco", dest='energyRes', default=True,
			action='store_false', help="Assume perfect energy reco.")
    parser.add_argument("--nfits", default=1, type=int,
			help="No. of pseudo-experiments to generate and fit.")
    args = parser.parse_args()

    mass = args.mass; distance = args.distance; nevts = args.nevents
    useTriggerEff = args.triggEff; useEnergyRes = args.energyRes

    spectrum = np.zeros((RESE-1)*REST).astype("float")

    createSpectrum(spectrum, mass, distance, nevts, useEnergyRes, useTriggerEff, noise)
    spectrum_max = np.max(spectrum)
    
    # arrays to store the events
    #eventEnergy = spectrum.intArray(int(events))
    #eventTime = spectrum.intArray(int(events))
    """
    eventEnergy = np.zeros(int(nevts)).astype("int32")
    eventTime = np.zeros(int(nevts)).astype("int32")
    """
    aa = datetime.now()
    for i in range(0,args.nfits):
        evts = sample_spectrum()
        eventEnergy = (np.array(evts).astype("int32"))[:,0]
        eventTime = (np.array(evts).astype("int32"))[:,1]
        # read in event
        #getEvent(eventEnergy, eventTime, mass, distance, nevts, i)
        # minimize
        x_min = minimize_scalar(llh, bounds=(0,mass+1.0), method='bounded', options={'disp': True, 'xatol':0.005})
        print("i=%4d | nfev: %2d | fit mass: %.5f | convergence: %5s |"
	      " message: %20s"%( i, x_min.nfev, x_min.x, x_min.success, x_min.message ))
        with open("DATA/masses_found_"+str(distance)+"Mpc_"+str(nevts)+"Events_"+str(mass)+"eV.txt", "a") as myfile:
            myfile.write(str(x_min.x) + '\n')

    bb = datetime.now()
    print "TIME", bb-aa

