import numpy as np
from ctypes import *

def get_funcs_to_expose():
	dll = CDLL('./spectrum.so', mode=RTLD_GLOBAL)
	func = dll.generateDist
	func2 = dll.getEvent
	func3 = dll.getLLH
	func4 = dll.createSpectrum
	func.argtypes = [c_double, c_double, c_double, POINTER(c_double), POINTER(c_double), c_bool]
	func2.argtypes = [POINTER(c_int), POINTER(c_int), c_double, c_double, c_double, c_int]
	func3.argtypes = [c_double, c_double, c_double, c_bool, c_bool, c_double, POINTER(c_int), POINTER(c_int)]
	func4.argtypes = [POINTER(c_double), c_double, c_double, c_double, c_bool, c_bool, c_double]
	return func, func2, func3, func4

__generateDist, __getEvent, __getLLH, __createSpectrum = get_funcs_to_expose()

def generateDist(mass, dist, nevts, distribution, triggerEffs, useEnergyRes):
	distribution_p = distribution.ctypes.data_as(POINTER(c_double))
	triggerEffs_p = triggerEffs.ctypes.data_as(POINTER(c_double))

	__generateDist(mass, dist, nevts, distribution_p, triggerEffs_p, useEnergyRes)

def getLLH(mass, distance, events, triggEff, energyRes, noise, eventTime, eventEnergy):
	eventTime_p = eventTime.ctypes.data_as(POINTER(c_int))
	eventEnergy_p = eventEnergy.ctypes.data_as(POINTER(c_int))
	return __getLLH(mass, distance, events, triggEff, energyRes, noise, eventTime_p, eventEnergy_p)

def getEvent(eventEnergy, eventTime, mass, distance, events, filenumber):
        eventTime_p = eventTime.ctypes.data_as(POINTER(c_int))
        eventEnergy_p = eventEnergy.ctypes.data_as(POINTER(c_int))
        __getEvent(eventEnergy_p, eventTime_p, mass, distance, events, filenumber)

def createSpectrum(spectrum, mass, distance, events, useEenergyRes, useTriggerEff, noise):
	spectrum_p = spectrum.ctypes.data_as(POINTER(c_double))
	__createSpectrum(spectrum_p, mass, distance, events, useEenergyRes, useTriggerEff, noise)

if __name__ == '__main__':
	mass = 1.0
	dist = 5.0
	nevts = 10.0
	distribution = np.zeros(600*1000).astype("float64")
	triggerEffs = np.ones(600).astype("float64")
	useEnergyRes = True

	generateDist(mass, dist, nevts, distribution, triggerEffs, useEnergyRes)
