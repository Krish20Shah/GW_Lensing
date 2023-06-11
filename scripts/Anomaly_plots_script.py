# mt = Total mass of BBH
# q = mass ratio
# l = inclination angle 

import numpy as np
import matplotlib.pyplot as plt
from pycbc.waveform import get_td_waveform,get_fd_waveform, fd_approximants
import pylab
import pycbc

mtotal = float(input("Enter the total Mass: \n"))
q = float(input("Enter the mass ratio : \n"))
l = float(input("Enter the value of inclination(in rad): "))

def h(mt,q,l):
    hp, hc = get_td_waveform(approximant="IMRPhenomXPHM",
                             mass1=mt / (1 + q), 
                             mass2=mt * q / (1 + q), 
                             delta_t=1.0/4096, 
                             f_lower=30.0, 
                             inclination=l)

    hpt = hp.to_frequencyseries(delta_f =0.01)

    shifted_hp_time , _ = get_td_waveform(approximant="IMRPhenomXPHM",
                             mass1=mt / (1 + q), 
                             mass2=mt * q / (1 + q), 
                             delta_t=1.0/4096, 
                             f_lower=30.0, 
                             coa_phase = -np.pi/4,
                             inclination=l)

    function = np.sqrt(1.0)*np.exp( -1j*np.pi/2)
    hpt *= function
    hp2t = hpt.to_timeseries(delta_t = hpt.delta_t)

    pylab.plot(hp2t.sample_times, hp2t,'-y', label=" Type II")
    pylab.plot(shifted_hp_time.sample_times, shifted_hp_time,':k', label=r'GR + $ \pi /4$')
    pylab.plot(hp.sample_times, hp,'-r', label='GR')

    pylab.xlim(-.25, .1)
    pylab.ylabel('Strain')
    pylab.xlabel('Time (s)')
    pylab.legend()
    pylab.savefig(f'Output_{mt}_{q}_{l}.jpg')



h(mtotal,q,l)