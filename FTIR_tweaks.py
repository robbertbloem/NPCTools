from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt


def import_FTIR_data(path, filename, delimiter = ","):

    paf = path + filename
    data = numpy.loadtxt(paf, delimiter = ",")
    
    # split array into frequency axis and intensity w
    w_axis = data[:,0]
    w = data[:,1]

    return w_axis, w, data
    
    
def tweak_FTIR(path, filename, n, delimiter = ","):

    
    
    # plot related
    fig = plt.figure(figsize = (12,8))
    n_ax = 2
    ax = [0] * n_ax
    ax_i = 0
    ax[ax_i] = fig.add_subplot(211)
    ax_i = 1
    ax[ax_i] = fig.add_subplot(212)
    
    w_axis, w, data = import_FTIR_data(path, filename, delimiter = ",")
    
    # split array into frequency axis and intensity w
    ax_i = 0
    ax[ax_i].plot(w_axis, w)
    
    # length of frequency data
    N_w = len(w) 
    
    # inverse FFT, zeropadded to twice the original length
    t = numpy.fft.ifft(w, n = 2 * N_w)
    
    # plot time domain, shift a copy of temp so that zero is in the middle
    ax_i = 1
    ax[ax_i].plot(numpy.real(numpy.fft.fftshift(t)))

    # make a new array for the time domain
    # subtract the average from the data to prevent a step
    t_new = numpy.zeros(len(t), dtype = numpy.cfloat)
    ave = numpy.average(t_new)
    t_new[:n] = t[:n] - ave 
    t_new[-n:] = t[-n:] - ave 

    # plot new time domain, shift a copy of temp so that zero is in the middle
    ax_i = 1
    ax[ax_i].plot(numpy.real(numpy.fft.fftshift(t_new)))

    # FFT, limit to first half
    w_new = numpy.fft.fft(t_new)[:N_w]

    # insert into array with wavelengths
    data[:,1] = numpy.real(w_new)

    # write to file
#     paf = path + "new_files/" + mess[m_i]["filename"]
#     numpy.savetxt(paf, data, delimiter = ",")   

    # plot related    
    ax_i = 0
    ax[ax_i].plot(w_axis, numpy.real(w_new))
    
#     a = numpy.real(w_new)[:]
#     x_axis = numpy.arange(N_w)
#     b = mess[m_i]["A"] * numpy.cos(mess[m_i]["B"] * numpy.pi * x_axis / (N_w) + mess[m_i]["D"]) + mess[m_i]["C"]
    
#     x = int(N_w/2)
#     
#     y = b[x]
#     b[-x:] = y
#     
#     print(b)
    
#     ax[ax_i].plot(freq_axis, b)
#     ax[ax_i].plot(freq_axis, numpy.real(w_new) - b)
    
#     ax[ax_i].plot(freq_axis, numpy.zeros(N_w), color = "black")
    
    ax_i = 0
#     if mess[m_i]["ylimw"] != (0,0):  
#         ax[ax_i].set_ylim(mess[m_i]["ylimw"])
    ax[ax_i].set_xlabel("Wavenumber (cm-1)")
    
    s = "%s, n = %i" % (filename, n)
    ax[ax_i].set_title(s)

    ax_i = 1
#     if mess[m_i]["ylimt"] != (0,0):  
#         ax[ax_i].set_ylim(mess[m_i]["ylimt"])
    ax[ax_i].set_xlabel("Time (indices)")
    
#     paf = path + "graphs/" + mess[m_i]["filename"][:-3] + "png"
#     
#     fig.savefig(paf)
    
    
    


plt.show()

if __name__ == "__main__": 
    pass