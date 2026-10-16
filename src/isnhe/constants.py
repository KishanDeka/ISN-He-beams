"""Physical constants and unit conventions.

Internal units
--------------
* distances: km
* velocities: km/s
* time: s
* ionization rate: s^-1
* phase-space density: cm^-3 (km/s)^-3
"""

GM_SUN = 1.32712440018e11      # km^3 s^-2, heliocentric gravitational parameter
AU_KM = 1.495978707e8          # km
KB = 1.380649e-23              # J/K
AMU = 1.66053906660e-27        # kg
EV = 1.602176634e-19           # J
M_HE = 4.002602 * AMU          # kg, helium atom mass
V_EARTH = 29.78                # km/s, mean orbital speed of Earth (circular orbit)
KM_TO_CM = 1.0e5               # cm/km

# Reference distance for normalising the ionization rate: beta(r) = beta0 * (R0 / r)^2
R0_KM = AU_KM
