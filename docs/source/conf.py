# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
import datetime
from importlib import import_module
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', )))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'H2Powerlaw'
copyright = '2024, Logan H. Jones'
author = 'Logan H. Jones'
release = '0.5.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


extensions = [
#'sphinx.ext.autodoc',
'numpydoc',
'sphinx.ext.autosummary',
]

#autosummary_generate = True
#autosummary_generate_overwrite = False
autoclass_content = 'both'
autodoc_typehints = "description"

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
