#!/usr/bin/python

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from obspy import read
from obspy.core.utcdatetime import UTCDateTime
from datetime import datetime
import pytz

net = 'AM' # network. usually AM "amateur"
sta = 'R0000' # station callsign. 5 uppercase alphanumeric characters beginning with R
ch = 'SHZ' # channel. should be either SHZ or EHZ depending on the model
loc = '00' # location. usually 00

now = UTCDateTime.now()
yr = now.year
day = now.strftime('%Y.%j')
skip = 60 * 5 # the number of seconds to skip back
short = now - (60*5) # datetime of analysis stary
tz = int(datetime.now(pytz.timezone('America/New_York')).strftime('%z'))/100
fmin = 0.1 # min frequency
fmax = 25 # max freq (should not exceed 25)
fminbp = 0.7 # lower bandpass limit
fmaxbp = 2 # upper bandpass limit

# output locations and filenames
outdir = '/opt/data/obs'
heli = os.path.join(outdir, day + '.png')
heliband = os.path.join(outdir, day + '-lowband.png')
spec = os.path.join(outdir, 'spec.png')
specband = os.path.join(outdir, 'spec-band.png')

mseedloc = '/opt/data/archive/' + str(yr) + '/' + net + '/' + sta + '/' + ch + '.D/' + net + '.' + sta + '.' + loc + '.' + ch + '.D.' + day

# read meta values from miniseed
st = read(mseedloc)
net = str(st[0].stats.network)
sta = str(st[0].stats.station)
loc = str(st[0].stats.location)
ch = str(st[0].stats.channel)
startt = str(st[0].stats.starttime)
sr = str(st[0].stats.sampling_rate)

sbp = st.copy() # copy for bandpass
sbp = sbp.filter('bandpass', freqmin=fminbp, freqmax=fmaxbp, zerophase=True)
spu = st.slice(starttime=short, endtime=now) # slice for main spectrogram
sps = sbp.slice(starttime=short, endtime=now) # slice for bandpass spectrogram

# make main and bandpass helicorders
st.plot(type='dayplot', title=net + '.' + sta + '.' + loc + '.' + ch + ' - ' + startt + ' - rate: ' + sr + 'Hz - range: 0-25Hz', vertical_scaling_range=8e3, outfile=heli, color=['k', 'r', 'b', 'g'], time_offset=tz)
sbp.plot(type='dayplot', title=net + '.' + sta + '.' + loc + '.' + ch + ' - ' + startt + ' - rate: ' + sr + 'Hz - bandpass: 0.7-2.0Hz', vertical_scaling_range=7e2, outfile=heliband, color=['k', 'r', 'b', 'g'], time_offset=tz)
del st
del sbp

# Plot the Spectrogram
# for reference-the old way of doing it:
# spec = sp.spectrogram(log=False, title='AM.RCB43 ' + str(sp[0].stats.starttime), outfile='/opt/data/obs/spec15.png', dbscale=False, cmap='viridis') #for reference

# filter (demean/constant)
sp = spu.detrend(type='constant')
ss = sps.detrend(type='constant')
del spu
del sps

## ---------------------------- ##
# make spectrogram figure 1
sfig1 = plt.figure(figsize=(16,6), dpi=100)
ax1 = sfig1.add_axes([0.068, 0.75, 0.85, 0.2]) #[left bottom width height]
ax2 = sfig1.add_axes([0.068, 0.1, 0.85, 0.6], sharex=ax1)
ax3 = sfig1.add_axes([0.931, 0.1, 0.03, 0.6])

# labels
startt = str(sp[0].stats.starttime)
endt = str(sp[0].stats.endtime)
sfig1.suptitle(net + '.' + sta + '.' + loc + '.' + ch + ' - ' + startt + '--' + endt + ' - samplerate: ' + sr + 'Hz - range: 0-25Hz')
ax1.set_ylabel('Counts')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Frequency [Hz]')
ax3.set_ylabel('Energy density [dimensionless]') # doesn't work

# make time vector
t = np.arange(sp[0].stats.npts) / sp[0].stats.sampling_rate

# plot waveform (top subfigure)
ax1.plot(t, sp[0].data, 'k')

# plot spectrogram (bottom subfigure)
sfig1 = sp[0].spectrogram(show=False, axes=ax2, log=False, dbscale=False, cmap='viridis')
mappable = ax2.images[0]
plt.colorbar(mappable=mappable, cax=ax3)

ax2.set_ylim(fmin, fmax)

plt.savefig(spec)

del sp

## ---------------------------- ##
# make spectrogram figure 2
sfig2 = plt.figure(figsize=(16,4), dpi=100)
ax1 = sfig2.add_axes([0.068, 0.600, 0.85, 0.3]) #[left bottom width height]
ax2 = sfig2.add_axes([0.068, 0.115, 0.85, 0.4], sharex=ax1)
ax3 = sfig2.add_axes([0.932, 0.115, 0.03, 0.4])

# labels
sfig2.suptitle(net + '.' + sta + '.' + loc + '.' + ch + ' - ' + startt + ' - samplerate: ' + sr + 'Hz - bandpass: 0.7-2.0Hz')
ax1.set_ylabel('Counts')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Frequency [Hz]')
ax3.set_ylabel('Energy density [dimensionless]') # doesn't work

# make time vector
t = np.arange(ss[0].stats.npts) / ss[0].stats.sampling_rate

# plot waveform (top subfigure)
ax1.plot(t, ss[0].data, 'k')

# plot spectrogram (bottom subfigure)
sfig2 = ss[0].spectrogram(show=False, axes=ax2, log=False, dbscale=False, cmap='viridis')
mappable = ax2.images[0]
plt.colorbar(mappable=mappable, cax=ax3)

ax2.set_ylim(fminbp, fmaxbp)

plt.savefig(os.path.join(specband)

del ss
