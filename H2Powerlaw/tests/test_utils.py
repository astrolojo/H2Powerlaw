## Data manipulation and fitting
import numpy as np
import astropy.units as u

## Visualization
from matplotlib.ticker import FormatStrFormatter

import pytest
from . import utils as utils
from . import H2Powerlaw

"""
All numerical tests use the observed fluxes for NGC 5033 and derived values from the tutorial notebook, with the exception of the H2 column and mass values, which use dummy values of PL slope and Tl.

As of 16 Sept 2026, tests should evaluate to the following numerical values:
test_obs_ratio : 2.25758134e+01 (Normalized empirical ratio of S0 to S1)
test_obs_ratio_uncert : 0.0808122 (Uncertainty of S1 in N ratio space)
test_nratio_model : 2.89885155 (Normalized theoretical ratio of S0 to S1 for n = 5., Tl = 100., j_norm = 1, Tu = 2000.)
test_do_fit : 4.72501703, 69.07622873 (Best-fit n and Tl in this Pythonic formulation of TS16)
test_total_column : ~12.7763 (Approximate N for n = 5., Tl = 100., omega = 1 arcsec)
test_calc_mass : ~20805375.6 (Approximate mass for  n = 5., Tl = 100., D = 14.8 Mpc)
"""

@pytest.fixture
def test_H2Model():
    p1 = H2Powerlaw.H2Model(flux =  1e-17 * np.array([3.66, 18.20, 6.35, 12.69]),
            flux_err = 1e-17 * np.array( [0.35, 1.04, 0.31, 1.91] ),
            j_obs = [0, 1, 2, 3], name = 'NGC 5033'
            )

    return p1
    

def test_is_jnorm_in():
    j_obs = np.array([0, 1, 2, 3])
    with pytest.raises(ValueError) as excinfo:
        utils._is_jnorm_in(j_obs, 9)
        
    assert "not in list" in str(excinfo.value)

def test_convert_to_cgs():
    x = 1 * u.W / u.m**2
    assert utils._convert_to_cgs(x) == 1000 ## Asserts that 1 W/m**2 = 1000 erg/s/cm**2

def test_obs_ratio(test_H2Model):
    p1 = test_H2Model

    x, j_x = p1.obs_ratio(p1.j_obs)

    assert (x[0] == pytest.approx(2.25758134e+01, rel = 1e-6))  & (j_x[0] == 'S0') 


def test_obs_ratio_uncert(test_H2Model):
    p1 = test_H2Model

    x, j_x = p1.obs_ratio_uncert(p1.j_obs)

    assert (x[1] == pytest.approx(0.0808122))  & (j_x[1] == 'S1') 

def test_nratio_model():
    j_obs = np.array([0, 1, 2, 3])
    x = H2Powerlaw.H2Model.nratio_model(j_obs, 5., 100, j_norm = 1, Tu = 2000.)

    assert (x[0] == pytest.approx(2.89885155)) & (1 + x[1] == pytest.approx(1.0))


def test_do_fit(test_H2Model):
    p1 = test_H2Model
    param, unc = p1.do_fit(p1.j_obs, overwrite = True, verbose = True)

    exp_slope, exp_Tl = 4.72501703, 69.07622873
    
    assert (param[0] == pytest.approx(exp_slope)) & (param[1] == pytest.approx(exp_Tl))


def test_total_column(test_H2Model):
    p1 = test_H2Model
    p1.omega = 1 * u.arcsec**2
    p1.distance = 14.8 * u.Mpc

    col = p1.total_column(p1.j_obs, 1, 5.0, 100.)

    assert np.log10(col) == pytest.approx(12.776318866405113) # For initial/dummy values of slope and Tl, you *should* get about 12.77.


def test_calc_mass(test_H2Model):
    p1 = test_H2Model
    p1.omega = 1 * u.arcsec**2
    p1.distance = 14.8 * u.Mpc

    m = p1.calc_mass(p1.j_obs, 1, 5., 100., p1.distance)

    assert (m.value == pytest.approx(20805375.6217)) & (m.unit == u.Msun)


"""
To-do: How to best test the plotting helpers?
"""
