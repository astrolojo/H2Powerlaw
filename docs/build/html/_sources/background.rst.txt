############
Background
############

Philosophy and Important Assumptions
=====================================

The H2 excitation model of `Togi & Smith (2016) <https://ui.adsabs.harvard.edu/abs/2016ApJ...830...18T/abstract>`_ (**TS16** hereafter and on other pages of these docs) is a simple framework for estimating molecular gas masses directly from H2 rotational lines in the mid-IR. Unlike methods which assume >2 gas components with discrete temperatures, the **TS16** model uses a power-law (PL) distribution of H2 temperatures in its excitation modeling.

The PL framework also allows the user to easily derive H2 gas masses between some temperature bounds. The upper bound ``Tu`` is set to 2000 K by default, since hotter gas contributes negligibly to the total mass. A justifiable choice of the lower bound ``Tl``, however, is critical to getting sensible gas masses - be sure to try a few values out!


A full account of the motivation, assumptions, and limitations of the **TS16** model is left to `the original paper <https://ui.adsabs.harvard.edu/abs/2016ApJ...830...18T/abstract>`_, but we outline a few highlights below:

* The model assumes ortho- and para-H2 (odd *J* and even *J*, respectively) are in local thermodynamic equilibrium (LTE). It further assumes a fixed ortho-to-para abundance ratio at any temperature. 

* The model **does not assume** any particular driver for the H2 excitation, although an excess of warm molecular gas (indicated through a shallow PL slope) could suggest that shocks and other collisional processes are important.

* The shape of the H2 temperature distribution at very cold temperatures is difficult to constrain. **TS16** showed that even with the low-energy S(0) transition to anchor the cooler gas, lower cutoff temperatures ``Tl`` < 80 K or so remain consistent with a target's best-fit PL. In other words, H2 gas below some temperature can contribute to the mass but not be excited enough to contribute to the line fluxes.

Differences from the IDL Code
==============================

While we have tried to replicate the original, IDL implementation of the **TS16** model with high fidelity, some small differences between the two were inevitable. In particular, we noted during development that the choice of minimization algorithm can slightly affect the best-fit power-law parameters and their uncertainties.

The original IDL code uses the Levenberg-Marquardt algorithm to do its fitting. However, the equivalent method in ``scipy.curve_fit`` doesn't allow for parameter bounds, which we find are necessary to prevent unphysical solutions. Instead, ``H2Powerlaw`` uses the Trust Region Reflective ('trf') algorithm in ``scipy.curve_fit``.

Putting the entire sample used in **TS16** into ``H2Powerlaw`` yields best-fit PL slopes that are generally consistent within uncertainties with the published values. However, we note that the parameter *uncertainties* in ``H2Powerlaw`` are systematically larger than those returned by the IDL code. This can propagate to your H2 mass estimates, especially if extrapolating to fairly cold (< 80 K) temperatures.