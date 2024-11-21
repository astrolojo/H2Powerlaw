## Data manipulation and fitting
import numpy as np
import astropy.units as u
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import astropy.constants as const
from operator import itemgetter

## Visualization
from matplotlib.ticker import FormatStrFormatter

## Self
#from H2Powerlaw import H2Model as model
    
def _evaluate_z(t):
    """
    Evaluate the partition function Z(T). To be used when calculating column density ratios.

    Parameters
    ----------
    t : float
        Temperature in Kelvin at which to evaluate Z(T).

    Returns
    -------
    zt : float
        Sum of the partition functions z_p and z_o (for para- and ortho-H2). Includes a term for the ground states of para-H2 (J = 0, g = 0) and ortho-H2 (J = 1, g = 9).

    """
    ## Partition function for para H<sub>2</sub>, even J_upper
    z_p = 1 + np.sum([g_lib[j] * np.exp(-1 * Ek_lib[j]/t) for j in j_lower_lib[::2]])
    
    ## Partition function for ortho H<sub>2</sub>, odd J_upper
    z_o = (9 * np.exp(-170./t))  + np.sum([g_lib[j] * np.exp(-1 * Ek_lib[j]/t) for j in j_lower_lib[1::2]])

    zt = z_p + z_o
        
    return zt

def _is_jnorm_in(j_obs, j_norm):
    """
    Checks that a given j_norm is allowable given the observed transitions.

    Parameters
    ----------
    j_norm : int
        List or numpy array of J_lower values associated with each line
    j_norm : int
        Value of J_lower to which all flux/column density ratios are normalized

    Returns
    -------
    None
    """

    if j_norm not in j_obs:
        raise ValueError("J_norm = " + str(j_norm) + " not in list of observed transitions, which include J = " + str(j_obs))

def _convert_to_cgs(model, arr):
    """
    Converts fluxes and flux uncertainties to cgs units (erg / s / cm^2) for computing H2 columns and masses. Note

    Parameters
    ----------
    arr : array-like
        List or array of flux/uncertainty values

    Returns
    -------
    arr_quantity : array-like
        Array of flux/uncertainty values after conversion to cgs units. Returns the value only.
    """

    arr_quantity = arr #* model.f_unit
    return arr_quantity.to(u.erg / u.s / u.cm**2).value



"""
The free-standing variables below define the energy levels and associated constants for 
the pure rotational levels of H2. For modeling purposes, the library currently only goes as high as the S(8) (J = 10 -> 8) transition but can be modified if needed. 

I define a similar library (which extends up to J_lower = 14), which currently only exists to
compute accurate partition functions. If/when I have the spoons to look up accurate wavelengths
and Einstein A coefficients, these could be consolidated into the other library described above.

"""

## Energy levels E_j / k (units Kelvin)
j_lower = np.arange(0,9,1)
j_upper = j_lower + 2
Ek = np.array( [510, 1015, 1681, 2503, 3473, 4585, 5828, 7196, 8677] )


## Corresponding wavelengths (microns) and Einstein A values (s^-1)
lam = np.array( [28.219, 17.035, 12.279, 9.665, 8.025, 6.91, 6.109, 5.511, 5.053] )
coeff_a = 1e-11 * np.array( [2.95, 47.6, 275., 980., 2640., 5880., 11400., 20000., 32400.] )

     
## Degeneracy values for para and ortho H2
## Access para H2 with g[::2] (even j) and ortho H2 with g[1::2] (odd j)
g = np.zeros(len(j_lower)) 

for j in j_upper:
    if j%2 == 0:
        g[j-2] += 2*j+1
    else:
        g[j-2] += 3*(2*j+1)
        
        
## A dictionary containing the above arrays. For easier access and indexing when some lines are absent
## (e.g., when S(0) isn't available, which may otherwise throw off the choice of normalization index j_norm)
linedict = dict()

for j in j_lower:
    tempdict = {"gu":g[j], "lam":lam[j], "A":coeff_a[j], "Eu":Ek[j]}
    
    linedict["S"+str(j)] = tempdict       


## Energy levels E_j / k (units Kelvin)
## Used only to accurately compute partition functions.
j_lower_lib = np.arange(0,15,1)
j_upper_lib = j_lower_lib + 2
Ek_lib = np.array( [510, 1015, 1681, 2503, 3473, 4585, 5828, 7196, 8677, 10263, 11940, 13703, 15549, 17458, 19402] )

g_lib = np.zeros(len(j_lower_lib)) ## Access para H2 with g[::2] (even j) and ortho H2 with g[1::2] (odd j)

for j in j_upper_lib:
    if j%2 == 0:
        g_lib[j-2] += 2*j+1
    else:
        g_lib[j-2] += 3*(2*j+1)
