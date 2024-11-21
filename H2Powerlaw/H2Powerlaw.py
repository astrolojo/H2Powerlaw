__all__ = [
    "linedict",
    "H2Model",
]

## Data manipulation and fitting
import numpy as np
import astropy.units as u
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import astropy.constants as const
from operator import itemgetter


## Visualization
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from IPython.display import display, Math

## Self
from . import utils as utils

linedict = utils.linedict

class H2Model:
    """The main API for the ``H2Powerlaw`` package, based on the Togi & Smith 2016 H2 excitation model.

    Users interact with data+model objects, which contain a set of empirical line fluxes, flux uncertainties, and the associated inputs and outputs for the power-law model. Not all parameters need to be set to perform the excitation fitting (e.g., distance) but will need to be set to compute a gas mass.

    Indexing of transition J values are done by mapping input arrays to a dictionary (imported from H2Powerlaw.utils) containing relevant line parameters. This should allow for arbitrary combinations of J values in the modeling without scrambling constants - e.g., if the user has J_lower = 1, 2, and 4, but not 3. If there's a smarter way to do this, please tell me.

    Parameters
    ----------
    flux : array-like
        List or numpy array of line fluxes
    flux_err : array-like
        List or numpy array of line flux uncertainties
    j_obs : array-like
        List or numpy array of J_lower values associated with each line
    j_norm : int, optional
        Value of J_lower for the transition to which other lines are normalized
    omega : astropy.units.quantity.Quantity, optional
        Beam size of observations, with units of solid angle
    distance : astropy.units.quantity.Quantity, optional
        Distance to object, with units of length/distance
    name : str, optional
        Name/label of object
    f_unit : astropy.units.core.CompositeUnit, optional
        Units of fluxes and uncertainties.

    """

    def __init__(self, flux, flux_err, j_obs, j_norm = 1, omega = None, distance = None, name = None, f_unit = u.W / u.m**2):

        self._flux = np.array(flux)
        self._flux_err = np.array(flux_err)
        self._funit = f_unit
        self._j_obs = np.array(j_obs)
        self._j_norm = j_norm
        self._omega = omega
        self._distance = distance
        self._Tu = 2000.
        self._Tl = 100.
        self._Tl_err = 0.
        self._slope = 5.
        self._slope_err = 0.
        self._name = str(name)

        # Ensure flux and flux_err parameters have sensible units
        if u.physical.get_physical_type(self._funit) != u.physical.energy_flux:
            raise u.UnitConversionError(str(self._funit)+" is not compatible with units of energy flux.")

        # Ensure distance parameter, if set, has sensible units
        if self._distance is not None:
            if type(self._distance) != u.quantity.Quantity:
                raise AttributeError("You may have forgotten to assign units to parameter 'distance.'")
            if u.physical.get_physical_type(self._distance) != u.physical.length:
                raise u.UnitConversionError(str(self._distance.unit)+" is not compatible with units of distance.")
        else:
            pass

        # Ensure beam size parameter, if set, has sensible units
        if self._omega is not None:
            if type(self._omega) != u.quantity.Quantity:
                raise AttributeError("You may have forgotten to assign units to parameter 'omega.'")
            if u.physical.get_physical_type(self._omega) != u.physical.solid_angle:
                raise u.UnitConversionError(str(self._omega.unit)+" is not compatible with units of solid angle.")
        else:
            pass


    """
    Define properties and setters, if needed
    """
    @property
    def flux(self):
        return self._flux
    @flux.getter
    def flux(self):
        return self._flux * self._funit
    @flux.setter
    def flux(self, array):
        self._flux = np.array(array)

    @property
    def flux_err(self):
        return self._flux_err
    @flux_err.getter
    def flux_err(self):
        return self._flux_err * self._funit
    @flux_err.setter
    def flux_err(self, array):
        self._flux_err = np.array(array)

    @property
    def f_unit(self):
        return self._funit
    @f_unit.setter
    def f_unit(self, units):
        if u.physical.get_physical_type(units) != u.physical.energy_flux:
            raise u.UnitConversionError(str(units)+" is not compatible with units of energy flux.")
        else:
            self._funit = units

    @property
    def j_norm(self):
        return self._j_norm
    @j_norm.setter
    def j_norm(self, value):
        self._is_jnorm_in(value)
        self._j_norm = value

    @property
    def Tu(self):
        return self._Tu
    @Tu.setter
    def Tu(self, value):
        self._Tu = value

    @property
    def Tl(self):
        return self._Tl
    @Tl.setter
    def Tl(self, value):
        self._Tl = value

    @property
    def slope(self):
        return self._slope
    @slope.setter
    def slope(self, value):
        self._slope = value

    @property
    def Tl_err(self):
        return self._Tl
    @Tl_err.setter
    def Tl_err(self, value):
        self._Tl_err = value

    @property
    def slope_err(self):
        return self._slope_err
    @slope_err.setter
    def slope_err(self, value):
        self._slope_err = value

    
    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, name):
        self._name = str(name)

    @property
    def j_obs(self):
        return self._j_obs
    @j_obs.setter
    def j_obs(self, array):
        self._j_obs = np.array(array) # Ensures j_obs is treated as a numpy array, even if input as a list
    
    @property
    def omega(self):
        return self._omega
    @omega.setter
    def omega(self, Quantity):
        if type(Quantity) != u.quantity.Quantity:
                raise AttributeError("You may have forgotten to assign units to parameter 'omega.'")
        if u.physical.get_physical_type(Quantity) != u.physical.solid_angle:
            raise u.UnitConversionError(str(Quantity.unit)+" is not compatible with units of solid angle.")
        else:
            self._omega = Quantity

    
    @property
    def distance(self):
        return self._distance
    @distance.setter
    def distance(self, Quantity):   
        if type(Quantity) != u.quantity.Quantity:
                raise AttributeError("You may have forgotten to assign units to parameter 'distance.'")
        if u.physical.get_physical_type(Quantity) != u.physical.length:
            raise u.UnitConversionError(str(Quantity.unit)+" is not compatible with units of length.")
        else:
            self._distance = Quantity

    def __str__(self):
        return f"H2 Excitation Model"

    
    def obs_ratio(self, j_obs):
        """
        Normalize observed fluxes (set when initializing the ``H2Model`` object) to some transition J_norm (also set at initialization). For flexibility, users specify J values of the transitions to use - this can be the entire set that was defined when initializing the model, or some subset. The chosen (sub)set MUST include the normalizing transition, default J_lower = 1.
        
        Parameters
        ----------
        j_obs : array-like
            List or numpy array of J_lower values associated with each line

        Returns
        -------
        fratios.values : dict values
            Array of flux ratios relative to J_norm
        fratios.keys : dict keys
            Array of J_lower values of chosen transitions. Used mainly for explicit tracking/indexing.
        """

        j_obs = np.array(j_obs)
        j_is_in = [j in j_obs for j in self.j_obs]
        flux_ = self.flux[j_is_in]
        flux_ = utils._convert_to_cgs(self, flux_)
        
        
        utils._is_jnorm_in(j_obs, self.j_norm)
        
                
        ## Initialize an empty dictionary for the column density ratios
        fratios = dict()
        
        ## Retrieve parameters relevant to transition j_norm
        g_norm, lam_norm, A_norm, Eu_norm = itemgetter('gu', 'lam', 'A', 'Eu')(linedict["S"+str(self.j_norm)])
        
    
        for j in utils.j_lower:
            if j not in j_obs:
                pass
            else:
                idx = np.where(j_obs == j)[0][0]
                idx_norm = np.where(j_obs == self.j_norm)[0][0]
                
                ## Retrieve parameters relevant to transition j
                g, lam, A, Eu = itemgetter('gu', 'lam', 'A', 'Eu')(linedict["S"+str(j)])
                
                ## Compute empirical column density ratios using TS16 eq. 11
                # fr_ = (g_norm / g) * ((self.flux[idx] * lam) / A) / \
                # ((self.flux[idx_norm] * lam_norm) / A_norm)
                fr_ = (g_norm / g) * ((flux_[idx] * lam) / A) / \
                ((flux_[idx_norm] * lam_norm) / A_norm)
                
                fratios["S"+str(j)] = fr_
    
        return np.array(list(fratios.values())), np.array(list(fratios.keys()))


    def obs_flux_uncert(self, j_obs):
        """
        Normalize flux uncertainties (set when initializing the ``H2Model`` object) to those for some transition j_norm (also set at initialization). For flexibility, users specify J values of the transitions to use - this can be the entire set that was defined when initializing the model, or some subset. The chosen (sub)set MUST include the normalizing transition, default J_lower = 1.
        
        Parameters
        ----------
        j_obs : array-like
            List or numpy array of J_lower values associated with each line

        Returns
        -------
        fr_uncert.values : dict.values
            Array of flux uncertainty ratios relative to J_norm
        fr_uncert.keys : dict.keys
            Array of J_lower values of chosen transitions. Used mainly for explicit tracking/indexing.
        """
        j_obs = np.array(j_obs)
        j_is_in = [j in j_obs for j in self.j_obs]
        flux_ = self.flux[j_is_in]
        flux_ = utils._convert_to_cgs(self, flux_)
        flux_err_ = self.flux_err[j_is_in]
        flux_err_ = utils._convert_to_cgs(self, flux_err_)

        utils._is_jnorm_in(j_obs, self.j_norm)
        
        ## Initialize an empty dictionary for the flux uncertainties
        fr_uncert = dict()
        
        ## Retrieve parameters relevant to transition j_norm
        g_norm, lam_norm, A_norm, Eu_norm = itemgetter('gu', 'lam', 'A', 'Eu')(linedict["S"+str(self.j_norm)])
    
        for j in utils.j_lower:
            if j not in j_obs:
                pass
            else:
                idx = np.where(j_obs == j)[0][0]
                idx_norm = np.where(j_obs == self.j_norm)[0][0]
                
                ## Retrieve parameters relevant to transition j
                g, lam, A, Eu = itemgetter('gu', 'lam', 'A', 'Eu')(linedict["S"+str(j)])
    
                ## Compute fractional flux uncertainty of transition j, normalized to j_norm
                ## Note that these are FLUX uncertainties, not column-density uncertainties
                # uncert_ratio_ = ( (self.flux_err[idx]/self.flux[idx])**2 + (self.flux_err[idx_norm]/self.flux[idx_norm])**2 )**0.5
                uncert_ratio_ = ( (flux_err_[idx]/flux_[idx])**2 + (flux_err_[idx_norm]/flux_[idx_norm])**2 )**0.5
                fr_uncert["S"+str(j)] = uncert_ratio_
    
        return np.array(list(fr_uncert.values())), np.array(list(fr_uncert.keys()))


    @classmethod
    def nratio_model(cls, j_obs, n, Tl, j_norm = 1, Tu = 2000.):
        """
        The theoretical column density ratios for transitions j_obs from the TS16 continuous temperature model (see their Equation 12). Can be used independently of a specific ``H2Model`` object to examine alternative fits.
        
        Parameters
        ----------
        j_obs : array-like
            List or numpy array of J_lower values associated with each line
        n : float
            Slope of the power-law temperature distribution. Fittable parameter
        Tl : float
            Lower temperature bound for the column density integral. Fittable parameter
        j_norm : int, optional
            Value of J_lower to which all flux/column density ratios are normalized
        Tu : float, optional
            Upper temperature bound for the column density integral

        Returns
        -------
        modelratios.values : dict.values
            Array of model ln(N) ratios relative to J_norm
        """

        j_obs = np.array(j_obs)
        ## Initialize an empty dictionary for the column density ratios
        modelratios = dict()
        
        ## Retrieve parameters relevant to transition j_norm
        g_norm, lam_norm, A_norm, Eu_norm = itemgetter('gu', 'lam', 'A', 'Eu')(linedict["S"+str(j_norm)])
        
        for j in utils.j_lower:
            if j not in j_obs:
                pass
                
            else:
                ## Retrieve parameters relevant to transition j
                g, lam, A, Eu = itemgetter('gu', 'lam', 'A', 'Eu')(linedict["S"+str(j)])
                
                ## Evaluate the numerator and denominator of TS16 equation 12
                numer = integrate.quad(lambda t: (1 / utils._evaluate_z(t)) * np.exp( -1 * Eu / t ) * t**(-1*n), Tl, Tu)
                denom = integrate.quad(lambda t: (1 / utils._evaluate_z(t)) * np.exp( -1 * Eu_norm / t ) * t**(-1*n), Tl, Tu)
                
                modelratios["S"+str(j)] = numer[0]/denom[0]
                
        return np.log( np.array(list(modelratios.values())) )


    def do_fit(self, j_obs, verbose = False, overwrite = False):
        """
        Fits the observed flux/column density ratios to the TS16 model for a given set of J_lower values. Uses scipy.curve_fit().
        
        Parameters
        ----------
        j_obs : array-like
            List or numpy array of J_lower values associated with each line
        verbose : bool, optional
            Provides numerical values for the observed and modeled column density ratios and fit results
        overwrite : bool, optional
            Sets whether to update the model values for the ``H2Model`` object

        Returns
        -------
        params : np.ndarray
            Array of best-fit slope and Tl, in that order
        params_uncert : np.ndarray
            Array of uncertainties on slope and Tl, in that order
        """

        utils._is_jnorm_in(j_obs, self.j_norm)
        
        nratio_obs, _ = self.obs_ratio(j_obs)
        nratio_obs_err, _ = self.obs_flux_uncert(j_obs)

        params, cov = curve_fit(self.nratio_model, j_obs, np.log(nratio_obs), 
                                sigma = nratio_obs_err, bounds = ([3., 20.], [7., 300.]))
        params_uncert = np.sqrt(np.diag(cov))

        print("Modeling column density ratios for  J_lower =",str(j_obs))
        if verbose == True:
            print("Nat. log of observed column density ratios, normalized to J =", self.j_norm, ":")
            print( np.log(nratio_obs) )
            print("Nat. log of model column density ratios, normalized to J =", self.j_norm, ":")
            print( self.nratio_model(j_obs, params[0], params[1]) )

        if overwrite == True:
            self.slope, self.Tl = params[0], params[1]
            self.slope_err, self.Tl_err = params_uncert[0], params_uncert[1]

        return params, params_uncert

    def total_column(self, j_obs, j_calc, n, Tl, Tu = 2000., log = False):
        """
        The total column density of H2(T > Tl), based on the flux of a given transition; see TS16 Eq 9. Depends strongly on the chosen value for beam size omega.
        
        Parameters
        ----------
        j_obs : array-like
            List or numpy array of J_lower values associated with each line
        j_calc : int
            Value of J_lower for the transition on which to base the column calculation
        n : float
            Slope of the power-law temperature distribution
        Tl : float
            Lower temperature bound for the column density integral
        Tu : float, optional
            Upper temperature bound for the column density integral
        log : bool, optional
            Sets whether to return N or log10(N)

        Returns
        -------
        column : float
            Column density of H2
        """

        if j_calc not in j_obs:
            raise Exception("Transition "+str(j_)+" not in list of observed lines.")
            
        else:
            j_obs = np.array(j_obs)
            j_is_in = [j in j_obs for j in self.j_obs]
            flux_ = self.flux[j_is_in]
            flux_ = utils._convert_to_cgs(self, flux_)
        
            idx = np.where(j_obs == j_calc)[0][0]
    
            ## Retrieve parameters relevant to transition j
            g, lam, A, Eu = itemgetter('gu', 'lam', 'A', 'Eu')(linedict["S"+str(j_calc)])
    
            ## Evaluate the definite integral in TS16 Equation 3
            n_integral = integrate.quad(lambda t: (g / utils._evaluate_z(t)) * np.exp( -1 * Eu / t ) * t**(-1*n), Tl, Tu)[0]
    
            ## Plug and chug values in TS16 Equation 9
            column = (4 * np.pi * flux_[idx] * lam * (Tl**(1 - n) - Tu**(1 - n))) / \
            (A * const.h.to("erg s").value * const.c.to("micron/s").value * self.omega.value * (n - 1) * n_integral)
    
            if log == True:
                return np.log10(column)
            else:
                return column

    def calc_mass(self, j_obs, j_calc, n, Tl, D, verbose = False, Tu = 2000., log = False):
        """
        The mass of H2(T > Tl), based on the flux of a given transition; see TS16 Eqs. 13 and 14.
        
        Parameters
        ----------
        j_obs : array-like
            List or numpy array of J_lower values associated with each line
        j_calc : int
            Value of J_lower for the transition on which to base the column calculation
        n : float
            Slope of the power-law temperature distribution
        Tl : float
            Lower temperature bound for the column density integral
        D : astropy.unit.quantity.Quantity
            Distance to target, with units specified
        verbose : bool, optional
            Prints mass information
        Tu : float, optional
            Upper temperature bound for the column density integral
        log : bool, optional
            Sets whether to return mass M or log(M)

        Returns
        -------
        mass_h2 : astropy.units.quantity.Quantity
            Mass of H2(T>Tl), in Msun
        """

        if type(D) != u.quantity.Quantity:
                raise AttributeError("You may have forgotten to assign units to distance parameter D.")
        if u.physical.get_physical_type(D) != u.physical.length:
            raise u.UnitConversionError(str(D.unit)+" is not compatible with units of length.")
        else:
            D = D.to(u.cm)
 
        
        total_column = self.total_column(j_obs, j_calc, n, Tl, Tu=Tu)
        total_number = total_column * self.omega.value * D.value**2
        mass_h2 = total_number * (3.32e-27 / const.M_sun.value)

        if log == True:
            m_ =  np.log10(mass_h2)
            if verbose == True:
                print("The total log(mass) of H2 between",round(Tl,1), "and", round(Tu,1),  "K is", round(m_, 2), "Msun")

        else:
            m_ = mass_h2
            if verbose == True:
                print("The total mass of H2 between",round(Tl,1), "and", round(Tu,1),  "K is", round(m_, 1), "Msun")

        return m_ * u.Msun

    def _plot_formatting(self, ax, xscale = 'linear'):
        """
        Small formatting details for excitation plots.
        
        Parameters
        ----------
        ax : pyplot Axis
            Axis on which to apply formatting
        xscale : str, optional ('linear' or 'log')
            Sets whether to plot the energy / x-axis in linear or log scale. 
        """

        ax.set_xlabel(r"E$_u$ / k (K)", weight='semibold')
        ax.set_ylabel(rf"ln(N$_u$ / g$_u$) - ln(N$_{self.j_norm}$ / g$_{self.j_norm}$)", weight='semibold')

        if xscale == "log":
            ax.set_xscale('log')
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        if self.name is not None:
            ax.set_title(self.name, weight='semibold')
            
        ax.legend(loc='best')

    def plot_obs_ratios(self, j_obs, obs_ratio,  ax, label = None, xscale = 'linear', **data_kwargs):
        """
        Plot the empirical column density ratios.
        
        Parameters
        ----------
        j_obs : array-like
            List or numpy array of J_lower values associated with each line
        obs_ratio : array-like
            List or numpy array of normalized column density ratios
        ax : pyplot Axis
            Axis on which to draw line
        label : str, optional
            Legend label; defaults to "Observed ratios"
        xscale : str, optional ('linear' or 'log')
            Sets whether to plot the energy / x-axis in linear or log scale. 
        **data_kwargs : dict
            Line and marker keyword arguments for matplotlib.pyplot.plot
        """

        E = [linedict['S'+str(j)]['Eu'] for j in j_obs]

        if label is None:
                label_ = "Observed ratios"
                obsline = ax.plot(E, np.log(obs_ratio),  
                        label = label_, **data_kwargs)
        else:
                obsline = ax.plot(E, np.log(obs_ratio),  
                        label = label, **data_kwargs)
        
        self._plot_formatting(ax, xscale)
        
        return obsline

    def plot_mod_ratios(self, j_obs, n, Tl, ax, label = None, **data_kwargs):
        """
        Plot the column density ratios from some power-law model.
        
        Parameters
        ----------
        j_obs : array-like
            List or numpy array of J_lower values associated with each line
        n : float
            Slope of the power-law temperature distribution
        Tl : float
            Lower temperature bound for the column density integral
        ax : pyplot Axis
            Axis on which to draw line
        label : str, optional
            Legend label; defaults to "Model ratios" followed by the given slope and temp.
        **data_kwargs : dict
            Line and marker keyword arguments for matplotlib.pyplot.plot
        """

        model_ratio = self.nratio_model(j_obs, n, Tl)
        E = [linedict['S'+str(j)]['Eu'] for j in j_obs]


        if label is None:
                label_ = "Model ratios (n= -"+str(round(n,2)) + ", Tl= " + str(round(Tl,1)) + ")"
                modline = ax.plot(E, model_ratio,  
                        label = label_, **data_kwargs)
        else:
                modline = ax.plot(E, model_ratio,  
                        label = label, **data_kwargs)
        
        self._plot_formatting(ax)
        
        return modline
