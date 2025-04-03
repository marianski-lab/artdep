# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import subprocess

for x in os.walk(os.path.abspath(os.path.join("..", "..", ".."))):
  sys.path.insert(0, x[0])

version_check = subprocess.run(['./version.sh'], stdout=subprocess.PIPE).stdout.decode('utf-8')
version_num = version_check.stdout

project = 'matmacore'
copyright = '2025, Eugene Chung, Ryan Kwok, Murat Yaman, Hillel Lerner, Mateusz Marianski'
author = 'Eugene Chung, Ryan Kwok, Murat Yaman, Hillel Lerner, Mateusz Marianski'
release = version_num

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              "sphinx.ext.githubpages",
              'sphinx.ext.todo',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

