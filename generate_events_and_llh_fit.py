from cppstuff_to_expose import getLLH, getEvent, createSpectrum
#import spectrum
from scipy.optimize import minimize
from scipy.optimize import brent
from scipy.optimize import minimize_scalar
from scipy.optimize import fmin
from scipy.stats import poisson
import numpy as np
from datetime import datetime
import argparse

def get_bin_centers(edges):
    return (edges[1:]+edges[:-1])/2.

# make binning configurable?
RESE = 200;
REST = 300;
EMAX = 60.0;
TMAX = 10.0;

def get_bin_edges():
    return [np.linspace(0, TMAX, REST+1), np.linspace(0, EMAX, RESE)]

def get_bin_centers(edgelist):
    return [(edges[1:]+edges[:-1])/2. for edges in edgelist]

bin_edges = get_bin_edges()
bin_centers = get_bin_centers(bin_edges)
T_bin_centers, E_bin_centers = bin_centers[0], bin_centers[1]

#triggEff = True
#energyRes = True
#mass = 1.0
#distance = 5.0
#events = 10
noise = 0.#pow(10,-3)*(10.0/REST)

# calculate llh for certain configuration
def getLLH(mass, distance, nevts, useEnergyRes, useTriggerRes, noise, pseudo_exp):
    hypo_spec = np.zeros((RESE-1)*REST).astype("float")
    createSpectrum(hypo_spec, mass, distance, nevts, useEnergyRes, useTriggerRes, noise)
    hypo_spec = hypo_spec.reshape(REST, RESE-1)
    if not np.alltrue(hypo_spec >= 0.0):
        raise ValueError(
            "Spectrum must have all bins >= 0.0! Spectrum generation bug?")
    val = -np.sum(poisson.logpmf(pseudo_exp, hypo_spec))
    return val

def llh(massi):
    #return getLLH( float(massi), distance, nevts, useTriggerEff, useEnergyRes, noise, eventTime, eventEnergy)
    return getLLH(massi, distance, nevts, useTriggerEff, useEnergyRes, noise, evts_histo)

def sample_spectrum():
    nevts_generated = 0
    evts = []
    #nsample = poisson.rvs(nevts)
    nsample = nevts
    while nevts_generated < nsample:
	# seed: make reproducible and non-redundant
        randE = np.random.randint(0, RESE-1)
        randT = np.random.randint(0, REST)
        rand_check = np.random.rand()*spectrum_max
        if spectrum[randT*(RESE-1)+randE] >= rand_check:
            #evts.append([randE, randT])
            evts.append([T_bin_centers[randT], E_bin_centers[randE]])
            nevts_generated += 1
    return np.array(evts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
		description="Simulate SN Spectrum and fit pseudo-experiments.")
    parser.add_argument("-m", "--mass", default=1.0, type=float,
			help="Neutrino mass to simulate")
    parser.add_argument("-d", "--distance", default=5.0, type=float,
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
    parser.add_argument("--asimov", default=False,
            action='store_true', help="Only fit the Asimov dataset.")
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
    if args.asimov:
        evts_histo = spectrum.reshape(REST, RESE-1)
        x_min = minimize_scalar(llh, bounds=(-mass-2.0, mass+2.0), method='bounded', options={'disp': True, 'xatol':0.005})
        print("Asimov | nfev: %2d | fit mass: %.5f | convergence: %5s |"
	      " message: %20s"%( x_min.nfev, x_min.x, x_min.success, x_min.message ))
    
    else:
        aa = datetime.now()
        for i in range(0,args.nfits):
            evts = sample_spectrum()
            np.savetxt("pseudo_exp_"+str(i)+"_"+str(distance)+"Mpc_"+str(nevts)+"Events_"+str(mass)+"eV.txt", evts)
            eventEnergy = evts[:,1]
            eventTime = evts[:,0]
            evts_histo, _, _ = np.histogram2d(eventTime, eventEnergy, bins=bin_edges)
            #evts_histo = poisson_sample().reshape(REST, RESE-1)
            # read in event
            #getEvent(eventEnergy, eventTime, mass, distance, nevts, i)
            # minimize
            x_min = minimize_scalar(llh, bounds=(0.0, 2.0), method='bounded', options={'disp': True, 'xatol':0.0000005})
            #x_min = minimize(llh, x0=[0.9], method='SLSQP', options={'disp': True, 'ftol': 1e-10})
            print("i=%4d | nfev: %2d | fit mass: %.5f | convergence: %5s |"
	          " message: %20s"%( i, x_min.nfev, x_min.x, x_min.success, x_min.message ))
            with open("DATA/masses_found_"+str(distance)+"Mpc_"+str(nevts)+"Events_"+str(mass)+"eV.txt", "a") as myfile:
                myfile.write(str(x_min.x) + '\n')

        bb = datetime.now()
        print "TIME", bb-aa

