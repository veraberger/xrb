import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.table import Table
from astropy.modeling import models
import glob
import scienceplots
import pdb

import pylag

from timing_funcs import *


# # ********* NICER *************
# # TODO read in evt file instead of a single gti lc from js
# # read in lc - currently a single gti
# cols = ['TIME(rel)', 'TIME(abs)', 'DEFAULT-RATE', 'NDETS', 'band0', 'band1', 'band2', 'band3', 'band4', 
#            'BGband5', 'BGband6', 'BGband7', 'BGband8', 'BGband9', 'HR1', 'HR2']
# # allGTIs = glob.glob('./*.lc')
# allGTIs = ['js_ni1200120131_0mpu7_silver_GTI19-v3-bands.lc']
# for f in allGTIs:
#     df = pd.read_csv(f, header=None, names=cols, skiprows=19, sep='\s+')
#     band0 = df['band0'] # 0.3-1 keV
#     band1 = df['band1'] #1-2
#     band2 = df['band2'] #2-4
#     band3 = df['band3'] #4-6
#     band4 = df['band4'] #6-12
#     t_x = df['TIME(abs)']
#     flux_x = band0+band1+band2+band3+band4

nicerfile = 'lc1200120131_0.001s.fits' # light curve output from xselect, 0.3-12 (13? oops?) kev

# # ********* HiPERCAM g band *************
paicefile = '/Users/vberger/xrb/paice_data_hcam_maxij1820_april172018/1820_Apr17_bary2_corrected_g.csv'
hlc = pd.read_csv(paicefile, header=0, skiprows=0, sep=',')
t_o = (hlc['Time_Bary']-hlc['Time_Bary'][0])*3600*24 # time in s, starting at t=0
flux_o = hlc['Target_g']


# pdb.set_trace()
# isolate the times that the nicer LC overlaps with the Hipercam times
# currently using one gti from js -> optical longer -> cut hcam lc
# flux_o = flux_o.loc[(t_o > t_x.iloc[0]) & (t_o < t_x.iloc[-1])]
# t_o = t_o.loc[(t_o > t_x.iloc[0]) & (t_o < t_x.iloc[-1])]

# # if x ray longer
# flux_x = flux_x.loc[(t_x > t_o.iloc[0]) & (t_x < t_o.iloc[-1])]
# t_x = t_x.loc[(t_x > t_o.iloc[0]) & (t_x < t_o.iloc[-1])]


print('read in data')

# # TIME LAG ANALYSIS W PYLAG
lc_o = pylag.LightCurve(t=t_o, r=flux_o) #  initializing pylag lc objects
lc_x = pylag.LightCurve(nicerfile) #interp_gaps=True)
print('made light curve objects')

# first hipercam obs is 58225.14557027361 MJD = 135401379.272 NICER time
# so subtract those times from both lcs and they will be aligned & hipercam lc starts at 0
lc_x.time -= 135401379.272

print('subtracted time')
# segments = lc_x.split_segments_time(segment_length=100) # split lc into segments of length 100 s
# print('segmented')


plt.plot(lc_x.time, lc_x.rate, label='x ray')
plt.plot(lc_o.time, lc_o.rate, label='optical')
plt.legend()
plt.show()

# REBIN LCs TO COMMON TIME RES [s]
"""
Mc lc and then put into equal length segments to take ffts of 
The length of those segments defines your lowest freq
Want enough data pts to take average and be in gaussian regime
"""
tbin = 0.003 # hcam time res [s]
# tbin = 0.125 # when working with js lc
binned_lc_x = lc_x.rebin(tbin=tbin)
print('rebinned x')
binned_lc_o = lc_o.rebin(tbin=tbin)


print('rebinned o, plotting....')
# pdb.set_trace()


# #sanity check
# plt.plot(t_o, flux_o)
# plt.plot(t_x, flux_x)
plt.plot(binned_lc_x.time, binned_lc_x.rate)
plt.plot(binned_lc_o.time, binned_lc_o.rate)
plt.show()

pdb.set_trace()


# #lc.time assumess, lc.rate assumes cts/s
# per = pylag.Periodogram(lc_x).bin(bins)
# pylag.Plot(per)

freq_o, ft_o= binned_lc_o.ft()
freq_x, ft_x = binned_lc_x.ft()
# freq_x, ft_x = lc_x.ft()
plt.plot(freq_x, ft_x)

plt.xscale('log')
plt.show()
# # sanity check
# plt.plot(freq_x, ft_x )
# plt.show()

cross_spec = pylag.CrossSpectrum(lc_o, lc_x)
# sanity check
frequencies = cross_spec.frequencies
cross_power = cross_spec.powers
plt.loglog(frequencies, abs(cross_power))
plt.xlabel("frequency (Hz)")
plt.ylabel("Cross Spectrum Power")
plt.show()



# # sanity check
# plt.plot(t_o, flux_o)
# plt.plot(t_x, flux_x)
# plt.show()

bins = pylag.LogBinning(1e-4, 1e-2, 20)

# lclist = get_lclist('*.lc')
# bins = pylag.LogBinning(1e-4, 1e-2, 20)
# stacked_per = pylag.StackedPeriodogram(lclist, bins)

# lc.zero_time() will reset the time axis to zero.
