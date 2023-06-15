import numpy as np
import matplotlib.pyplot as plt
from pycbc.waveform import get_td_waveform,get_fd_waveform, fd_approximants
import pylab
import pycbc



q = 0.3
mt = pycbc.conversions.mass2_from_mchirp_q(26.1,q)*(1+q)
l = [0,np.pi/3]
spin1 = [0.5,0,0]
spin2 = [0,0,0.5]

for i in range(len(l)):
    hpa, hca = get_td_waveform(approximant="IMRPhenomXPHM",
                             mass1=mt / (1 + q), 
                             mass2=mt * q / (1 + q), 
                             delta_t=1.0/4096, 
                             f_lower=30.0, 
                             inclination=l[i],
                             spin1x = spin1[0],spin1y = spin1[1],spin1z = spin1[2],
                             spin2x = spin2[0],spin2y = spin2[1],spin2z = spin2[2],
                             distance = 410)
    ht = hpa.to_frequencyseries()
    function = np.sqrt(1.0)*np.exp( -1j*np.pi/2)
    ht *= function
    h2t = ht.to_timeseries(delta_t = ht.delta_t)



    shifted_hpa, shifted_hca = get_td_waveform(approximant="IMRPhenomXPHM",
                             mass1=mt / (1 + q), 
                             mass2=mt * q / (1 + q), 
                             delta_t=1.0/4096, 
                             f_lower=30.0, 
                             inclination=l[i],
                             coa_phase = -np.pi/4,
                             delta_f = 0.01,
                             spin1x = spin1[0],spin1y = spin1[1],spin1z = spin1[2],
                             spin2x = spin2[0],spin2y = spin2[1],spin2z = spin2[2],
                             distance = 410)




    pylab.subplot(2*len(l),1,i+1)
    #pylab.plot(h.sample_times, h,'-r', label=" GR -a")
    pylab.plot(shifted_hpa.sample_times, shifted_hpa,'y', label=r'GR + $ \pi /4$' 'inclination = %1.3f'%l[i])
    pylab.plot(h2t.sample_times, h2t,':b', label='Type-II inclination = %1.3f'%l[i])
    #pylab.plot(hpm.sample_times, hpm,'-y', label='GR-m')
    #pylab.plot(hpn.sample_times, hpn,'-g', label='GR-n')
    pylab.xlim(-.25, .1)
    pylab.ylabel('Strain')
    pylab.xlabel('Time (s)')
    pylab.legend(loc = "lower right")

    h_residual = h2t-shifted_hpa


    pylab.subplot(2*len(l),1,i+1+len(l))
    pylab.plot(h_residual.sample_times, h_residual,'-k', label='Residual inclination = %1.3f'%l[i])
    pylab.xlim(-.25, .1)
    pylab.ylabel('Strain')
    pylab.xlabel('Time (s)')
    pylab.legend(loc = "lower right")    
plt.show()





