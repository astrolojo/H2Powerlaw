############
Quickstart
############

Overview and Dependencies
=========================

``H2Powerlaw`` model objects contain both a set of observed data and the parameters needed to generate a theoretical power-law (PL) excitation curve. For the example code below, we use H2 line fluxes for the galaxy NGC 5033 (given in Table 2 of Togi & Smith 2016, which in turn are obtained from the `SINGS <https://ui.adsabs.harvard.edu/abs/2003PASP..115..928K/abstract>`_ program). 

This example is also provided as a notebook, which can be downloaded from the `git repo <https://github.com/astrolojo/SESAMME/issues>`_. 

``H2Powerlaw`` relies on the following packages for its functionality. If you do not have these packages installed, you can install them using pip or conda.

* **numpy** and **scipy** for handling arrays, including fitting
* **astropy** for unit manipulation
* **matplotlib** for data visualization 


Using the ``H2Model`` Class
===========================

We start by importing the ``H2Model`` class, which contains all the important methods for doing your excitation fitting, and initializing our object. At minimum, every ``H2Model`` requires fluxes, flux uncertainties, and *J* values for the relevant transitions to perform an excitation fit, plus units of the flux measurements.

.. code-block:: python

 from H2Powerlaw import H2Model
 p1 = H2Model(flux =  1e-17 * np.array([3.66, 18.20, 06.35, 12.69]),
            flux_err = 1e-17 * np.array( [0.35, 1.04, 0.31,1.91] ),
            j_obs = [0,1,2,3],
            f_unit = u.W / u.m**2
            )

By default, fluxes are assumed to be in SI units of W/m^2, but any compatible units are converted to cgs (erg/s/cm^2) under the hood.

The other key parameter is ``j_norm``, which denotes which line to use as our normalizing transition when computing flux/column density ratios. By default ``j_norm = 1``, consistent with **TS16**.

Creating models and fitting
----------------------------

Column density ratios for a given PL model can be generated with the ``nratio_model()`` class method, meaning they can be created independently of any observed dataset or ``H2Model`` object. Given the same input arguments, both of the calls below will produce the same output.

.. code-block:: python

 slope, Tl = 5.1, 84.0 # Example PL index slope and lower temperature bound Tl
 
 p1.nratio_model(p1.j_obs, slope, Tl)
 H2Model.nratio_model([0 , 1, 2, 3], slope, Tl)


To fit a PL model to real data, simply call on ``H2Model.do_fit()``:

.. code-block:: python

 param, unc = p1.do_fit(p1.j_obs)

This returns the parameters (``param = (slope, Tl)``) and parameter uncertainties (``unc = (slope_err, Tl_err)``) for the best-fit PL model normalized to some ``j_norm``. You may also print the natural log of the observed and model column density ratios by setting ``verbose = True``.

Best-fit values and uncertainties can be stored through the properties ``p1.slope``, ``p1.slope_err``, etc. by setting ``overwrite = True``.


Note on optional parameters
-----------------------------
Other keyword arguments for ``H2Model`` describe the target object and the observations themselves. While not necessary to perform an excitation fit, they're critical to extracting H2 mass measurements in the **TS16** framework.

.. code-block:: python

 p1.name = "NGC 5033" # Target name
 p1.distance = 14.8 * u.Mpc # Distance to target, in Mpc
 p1.omega = 1 * u.arcsec**2 # Beam size of observations


Computing gas masses
----------------------

The PL model can also yield an estimate of the column density of H2 gas warmer than ``Tl`` (and less than some upper bound ``Tu``, default set to 2000 K):

.. code-block:: python

 col = p1.total_column(p1.j_obs, 1, p1.slope, p1.Tl)

Similarly, you can translate this column density to a mass at T > ``Tl``:

.. code-block:: python

 m = p1.calc_mass(p1.j_obs, 1, p1.slope, p1.Tl, p1.distance)


You may also explicitly print the H2 mass (in Msun) and temperature bounds by setting ``verbose = True``. 
