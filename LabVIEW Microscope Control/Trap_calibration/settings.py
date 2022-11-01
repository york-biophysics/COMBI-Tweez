from numpy import pi

# http://biopt.ub.edu/force-detection/brownian

kT = 1.381e-23 * 293.0
Fs = 100e3
r = 1.5  # bead radius in microns
gamma = 3 * pi * 8.9e-4 * 2*r*1e-6  # drag coefficient= 3pi*eta*d
D = kT/gamma
h = 1.0     # height above surface in microns

RoH = r/(r+h)

lat_correction = 1 - (9/16 * RoH) - (1/8 * RoH**3) - \
    (45/256 * RoH**4) - (1/16 * RoH**5)

axl_correction = 1 - (9/8 * RoH) + (1/2 * RoH**3) - (57/100 * RoH**4) - \
    (1/5 * RoH**5) - (7/200 * RoH**11) - (1/25 * RoH**12)
