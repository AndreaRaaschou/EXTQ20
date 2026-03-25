"""This file defines variables that can be used to calculate stroke work.
It is part of a hand-in task in the course EXTQ20

Copyright Einar Heiberg. For information or questions, please contact
Einar Heiberg, einar.heiberg@med.lu.se.

Data extracted from a healthy volunteer when not stated otherwise.
Software used for data extraction was Segment (http://segment.heiberg.se),
freely available for research and educational use.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from polyarea import polyarea

# General data
rho = 1050  # Density of blood [g/dm^3, or kg/m^3]

heartrate = 60  # Beats per minute. Actually true heart rate was 56, but use 60 to simplify things.

meanaorticvalvearea = 6.46  # cm^2

mmHg2Nperm2 = 133.3  # 1 mmHg is equal to 133.3 N/m^2, good to know conversion.

# Later interpolate to equidistant time interval, use 100 timesteps
nframes = 100

# Volume curve of the left ventricle volume. Data taken from
# healthy volunteer.
lvv = np.array([
    172.3350, 168.8071, 164.2898, 156.0784,
    141.7951, 123.5322, 104.4061, 89.4800,
    79.1059, 73.2753, 69.4553, 67.8909,
    71.1016, 81.5572, 98.2768, 117.2011,
    131.4172, 139.2867, 142.2792, 143.8322,
    145.0454, 145.9041, 146.6225, 147.4332,
    148.0553, 148.4486, 149.3266, 153.1151,
    161.8366, 170.4185
])

# Mean velocity in the proximal aorta. Velocity
# in cm/s.
vel = np.array([
    -1.358230, 1.016271, 10.293329, 30.930073,
    56.763535, 66.770241, 71.462341, 68.551804,
    63.509338, 56.910000, 48.758038, 38.965691,
    27.829584, 12.282402, -3.711286, -5.208059,
    -3.562259, -1.527832, 1.650186, 2.686651,
    3.293284, 3.190372, 3.405491, 3.260243,
    3.535447, 3.332693, 2.627509, 2.042885,
    1.245768, 1.358219, 2.086271, 2.482973,
    1.676623, 1.351324, 0.860088, 1.356620,
    1.899150, 3.266571, 3.034079, 3.530154
])

# Mean flow through the aortic valve in ml/s.
meanflow = np.array([
    -8.066255, 6.676226, 65.507202, 195.026154,
    371.234009, 438.636108, 480.987061, 467.428619,
    437.702362, 395.558197, 330.317535, 267.977692,
    189.351242, 80.687225, -23.999735, -33.908012,
    -23.349440, -9.992037, 10.816427, 17.491928,
    21.344938, 20.584368, 21.972319, 21.322018,
    22.707109, 21.258274, 16.605993, 12.881185,
    7.855052, 8.464514, 13.154749, 15.801749,
    10.694687, 8.540438, 5.208784, 8.195945,
    12.587546, 21.267601, 18.864103, 21.275459
])

# Mean left ventricle pressure (inside the ventricle).
# Define it as a sort of spline curve.
S = np.array([
    [0.00, 7, 85],
    [0.02, 10, 83],
    [0.05, 60, 80],
    [0.13, 115, 114],
    [0.20, 120, 120],
    [0.25, 115, 116],
    [0.30, 109, 111],
    [0.35, 100, 104],
    [0.38, 80, 107],
    [0.40, 40, 103],
    [0.43, 20, 102],
    [0.50, 8, 100],
    [0.60, 3, 90],
    [0.70, 2, 85],
    [0.95, 3, 84],
    [1.00, 7, 85]
])

# Left ventricle pressure
lvp_interp = PchipInterpolator(S[:, 0], S[:, 1])
lvp = lvp_interp(np.linspace(0, 1, nframes))  # result in mmHg

# Left ventricular aortic pressure
lvap_interp = PchipInterpolator(S[:, 0], S[:, 2])
lvap = lvap_interp(np.linspace(0, 1, nframes))  # result in mmHg

# Define sampling times
t = np.linspace(0, 1, nframes)  # s
deltat = t[1]  # It is 1/100 s between each timestep.

# Interpolate all signals to nframes.
lvv_interp = PchipInterpolator(np.arange(len(lvv)), lvv)
lvv = lvv_interp(np.linspace(0, len(lvv) - 1, nframes))

vel_interp = PchipInterpolator(np.arange(len(vel)), vel)
vel = vel_interp(np.linspace(0, len(vel) - 1, nframes))

meanflow_interp = PchipInterpolator(np.arange(len(meanflow)), meanflow)
meanflow = meanflow_interp(np.linspace(0, len(meanflow) - 1, nframes))

# Task 2
# Calculate the derivative of the lvv
x_values  = np.arange(0, 100)
dlvvdt = np.array([(lvv[x] - lvv[x-1]) / deltat for x in x_values])

# Task 1
# Make plots of the given data in strokework
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, num=10, figsize=(10, 8))

# left ventricle volume with EDV and ESV
ax1.plot(t, lvv, 'k.-', label = 'LV Volume [ml]')
ax1.plot(t[39], lvv[39], 'o', label = 'End systolic volume (ESV)')
ax1.plot(t[99], lvv[99], 'o', label = 'End diastolic volume (EDV)')
ax1.grid(True)
ax1.legend()
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Volume [ml]')
ax1.set_title('LV volume [ml]')

# combined: pressure curve in the aorta and the left ventricle (lvap and lvp)
ax2.plot(t, lvp, 'b.-', label='LV pressure [mmHg]')
ax2.plot(t, lvap, 'r.-', label='Aortic Pressure [mmHg]')
ax2.grid(True)
ax2.set_title('Pressure curves')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Pressure [mmHg]')
ax2.legend(loc='upper right')

# combined: aortic flow (meanflow) and dLVV/dt
ax3.plot(t, meanflow, 'k.-', label = 'Aortic flow curve')
ax3.plot(t, dlvvdt, 'r.-', label = 'dLVV/dt')
ax3.plot(t[dlvvdt.argmax()], np.max(dlvvdt), 'o', label = 'Peak filling rate time (PFRT)')
ax3.plot(t[dlvvdt.argmin()], np.min(dlvvdt), 'go', label = 'Peak ejection rate time (PERT)')
ax3.grid(True)
ax3.set_title('Aortic mean flow curve [ml/s]')
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Flow [ml/s]')
ax3.legend(loc='upper right')

# velocity curve (vel)
ax4.plot(t, vel, 'k.-')
ax4.grid(True)
ax4.set_title('Velocity curve of blood in the aorta')
ax4.set_xlabel('Time [s]')
ax4.set_ylabel('Pressure [cm]')


# user interaction to get aortic valve opening time and pressure
#print("Time and pressure of the valve opening and closing:")
#x = plt.ginput(2)
#print(x)

plt.tight_layout()
plt.show()

# Calculate SV, CO and EF
SV = lvv[99] - lvv[39]
CO = SV*60/1000 # considering a heart rate of 60 beats per min
EF = SV / lvv[99]
print(f'Stroke volume (SV) = {SV} ml')
print(f'Cardiac output (CO) = {CO} l/min')
print(f'Ejection fraction (EF) = {EF}')

# Task 4: Plotting the pressure volume loop
plt.plot(lvv, lvp, 'k.-', label = 'PV loop')
plt.title( 'Pressure volume (PV) loop')
plt.xlabel('Volume (ml)')
plt.ylabel( 'Pressure (mmHg)')
plt.plot(lvv[39], lvp[39], 'bo', label = 'ESV')
plt.plot(lvv[99], lvp[99], 'ro', label = 'EDV')
plt.plot(lvv[dlvvdt.argmax()], lvp[dlvvdt.argmax()],'go', label = 'PFRT')
plt.plot(lvv[dlvvdt.argmin()], lvp[dlvvdt.argmin()], 'yo', label = 'PERT')
plt.legend(loc = 'upper right')

# Task 4: compute the area of the polygon - stroke work (SW)
SW_1 = polyarea(lvv, lvp) * 1.33322e-4
SW_2 = abs(np.trapezoid(lvp* dlvvdt, t)) * 1.33322e-4
SW_3 = SV * np.mean(lvap) * 1.33322e-4
print(f'Stroke work (SW_1) = {SW_1:.3f} J')
print(f'Stroke work (SW_2) = {SW_2:.3f} J')
print(f'Stroke work (SW_3) = {SW_3:.3f} J')

# Task 5: Compute the KE
KE = 0.5 * SV*1.05 * 10**(-3) * (10**(-2) * np.mean(meanflow[0:40])/meanaorticvalvearea)**2
print(f'Kinetic energy (KE) = {KE:.5f} J')

# Task 6: Change in SW and KE with increased activity
SW_inc = SV * 1.3 * np.mean(lvap)*1.5 * 1.33322e-4
print(f'Stroke work during exercise (SW_inc) = {SW_inc:.3f} J')

KE_inc = 0.5 * SV*1.3*1.05 * 10**(-3) * (10**(-2) * np.mean(meanflow[0:40])*2.2/meanaorticvalvearea)**2
print(f'Kinetic energy during exercise (KE_inc) = {KE_inc:.5f} J')

#print(np.mean(meanflow[0:40]))

print(KE_inc / KE)
print(SW_inc / SW_3)














