# TWM-Semiconductor-Laser

## Description of the Model
This repository contains a Travelling-Wave-Model (TWM) simulation for a simple semiconductor laser. The laser consists of a gain section enclosed between two mirrors with reflectivities _r1_ and _r2_. The simulation follows the equation model proposed by **A. Vladimirov et al.** in the paper **__"Numerical Study of Dynamical Regimes in a Monolithic Passively Mode-Locked Semiconductor Laser", IEEE Journal of Quantum Electronics, Vol. 45, No. 5, May 2009__**.

Key differences between the model proposed in the paper and the Python simulation are the absence of a saturable absorber and a passive section, as well as the absence of a thin spectral filter as described by _Eq.(4)_. Instead of the thin spectral filter, the same boundary condition, described by _Eq.(3)_, has been applied to both sides.

## Future Outlook
The simulation needs to be refined and updated to include all components as described in Vladimirov's paper. 
