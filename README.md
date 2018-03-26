# stellar_calculations
This code can be used to make approximate stellar structure calculations.

   General Description:
   ====================
        StatStar computes a ZAMS model using a number of simplifying assumptions about the physics.  The code
        is not designed to produce precise, research-quality models of ZAMS stars; rather, it is meant to 
        illustrate many of the fundamental physical ideas discussed in
 
            "An Introduction to Modern Astrophysics"
            Bradley W. Carroll and Dale A. Ostlie
            Second Edition, Addison Wesley,   Copyright 2007.
 
        StatStar performs an inward integration of the stellar structure equations, given values for the
        star's mass, luminosity, effective temperature, and composition.
 
        The simplifying assumptions made here include:
            (a) A constant composition throughout (characteristic of a ZAMS star).
            (b) The use of the ideal gas law throughout.
            (c) The gas is assumed to be completely ionized throughout.
            (d) Radiation pressure is incorporated.
            (e) Convection, when present, is taken to be purely adiabatic.
            (f) Opacity is computed by using approximation formulae:
                1.  Bound-free processes via Kramer's formula.
                2.  Free-free processes via Kramer's formula.
                3.  Electron scattering via Thomson formula.
                4.  H- ion via fitting function.
            (g) Surface is assumed to have P = 0, T = 0, rho = 0.
            (h) Outermost (optically thin) zone is assumed to be radiative.
            (i) No attempt is made to satisfy the Eddington approximation by
                adjusting the outside boundary condition.
