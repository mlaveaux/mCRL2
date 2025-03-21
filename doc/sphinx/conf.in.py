# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Imports -----------------------------------------------------------------

import datetime
import os
import re
import sys

# used to parse mCRL2 version information from CMake file
from pathlib import Path

_CMAKE_SOURCE_DIR = '@CMAKE_SOURCE_DIR@'
_MCRL2_TOOL_PATH = '@MCRL2_TOOL_PATH@'
_MCRL2_TOOLS = '@MCRL2_TOOLS@'
_SPHINX_BUILD_TEMP_DIR = '@SPHINX_BUILD_TEMP_DIR@'
_SPHINX_BUILD_OUT_DIR = '@SPHINX_BUILD_OUT_DIR@'

os.environ['CMAKE_SOURCE_DIR'] = _CMAKE_SOURCE_DIR

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
sys.path.insert(0, (Path(__file__).parent / '_extensions').as_posix())

# -- Project information -----------------------------------------------------
project = 'mCRL2'
author = 'Technische Universiteit Eindhoven'
release = '' # populated below
version = ''

# run CMake on the version file to obtain the current version of mCRL2
from manual import call
olddir = os.getcwd()
try:
    os.chdir(_CMAKE_SOURCE_DIR)
    out = call('CMake', ['cmake', '-P',
        f'{_CMAKE_SOURCE_DIR}/cmake/MCRL2Version.cmake']).decode('utf-8')
    for line in iter(out.splitlines()):
        matches = re.findall(r'MCRL2_MAJOR_VERSION ([\S]+)', line)
        if matches:
            release = matches[0]
        matches = re.findall(r'MCRL2_MINOR_VERSION ([\S]+)', line)
        if matches:
            version = matches[0]
finally:
    os.chdir(olddir)
version = f'{release}.{version}'

# update copyright from current version
copyright = f'2011-{datetime.datetime.now().year}, {author}'


# -- General configuration ---------------------------------------------------

# The master document is in the user_manual directory.
master_doc = 'index'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'mcrl2_dparser',
    'mcrl2_pygment',
    'sphinx.ext.ifconfig',
    'sphinx.ext.imgmath',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme',
]

imgmath_image_format = 'svg'
imgmath_font_size = 14
imgmath_latex_preamble = '\\usepackage{newtxsf}'

myst_enable_extensions = [
    "dollarmath"
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['.git']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_style = 'css/mcrl2.css'
html_theme_options = {
    'collapse_navigation': False,
    'logo_only': True,
    'style_external_links': False
}
html_logo = '_static/mcrl2.png'
html_sidebars = {
   '**': ['localtoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html'],
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_extra_path = ['_doxygen/output/']
html_js_files = []

# Tweaking how the "last updated" is displayed
today_fmt = '%d-%m-%Y'

# Suppress warnings about unreferenced citations (there are many on the publications page).
suppress_warnings = ['ref.citation']

if tags.has('build_doxygen'):
    extensions.append('sphinxcontrib.doxylink')

doxylink = {
    'mcrl2' : (str(Path(__file__).parent / '_doxygen/mcrl2.tag'), 'doxygen/')
}
    
# -- App setup - executed before the build process starts --------------------
def setup(app):
    import manual
    import pdflatex

    os.makedirs(_SPHINX_BUILD_TEMP_DIR, mode = 0o755, exist_ok = True)
    os.makedirs(_SPHINX_BUILD_OUT_DIR, mode = 0o755, exist_ok = True)

    olddir = os.getcwd()
    try:
        os.chdir(_SPHINX_BUILD_TEMP_DIR)

        if tags.has('build_pdflatex'):
            pdflatex.generate_pdfs()

        if tags.has('build_manual'):
            manual.generate_rst(_CMAKE_SOURCE_DIR, _MCRL2_TOOL_PATH, _MCRL2_TOOLS)
    finally:
        os.chdir(olddir)
