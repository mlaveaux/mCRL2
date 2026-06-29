# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Imports -----------------------------------------------------------------

import datetime
import os
import re
import shutil
import sys

# used to parse mCRL2 version information from CMake file
from pathlib import Path

_CMAKE_SOURCE_DIR = '@CMAKE_SOURCE_DIR@'
_SPHINX_BUILD_TEMP_DIR = '@SPHINX_BUILD_TEMP_DIR@'
_SPHINX_BUILD_OUT_DIR = '@SPHINX_BUILD_OUT_DIR@'

os.environ['CMAKE_SOURCE_DIR'] = _CMAKE_SOURCE_DIR
os.environ['MCRL2_TOOLS'] = '@MCRL2_TOOLS@'
os.environ['MCRL2_TOOL_PATH'] = '@MCRL2_TOOL_PATH@'

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
sys.path.insert(0, (Path(__file__).parent / '_extensions').as_posix())

# -- Project information -----------------------------------------------------
project = 'mCRL2'
author = 'Technische Universiteit Eindhoven'
release = '@MCRL2_MAJOR_VERSION@'
version = '@MCRL2_MAJOR_VERSION@.@MCRL2_MINOR_VERSION@'

# update copyright from current version
copyright = f'2011-{datetime.datetime.now().year}, {author}'


# -- General configuration ---------------------------------------------------

# The master document is in the user_manual directory.
master_doc = 'index'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'mcrl2_admonitions',
    'mcrl2_dparser',
    'mcrl2_pygment',
    'sphinx.ext.ifconfig',
    'sphinx.ext.imgmath',
    'mcrl2_imgmath',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme',
]

imgmath_font_size = 14
imgmath_image_format = 'svg'
imgmath_latex_preamble = '\\usepackage{mcrl2_package}'
imgmath_use_preview = True

myst_enable_extensions = [
    "dollarmath",
    "colon_fence"
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
html_js_files = []

# Tweaking how the "last updated" is displayed
today_fmt = '%d-%m-%Y'

# Enable numbered figures and tables so :numref: works
numfig = True

# Suppress warnings about unreferenced citations (there are many on the publications page).
suppress_warnings = ['ref.citation']

if tags.has('build_doxygen'):
    extensions.append('sphinxcontrib.doxylink')
    html_extra_path = ['@DOXYGEN_OUTPUT_PATH@']
    doxylink = {
        'mcrl2' : ('@DOXYGEN_TAG_PATH@', 'doxygen/')
    }

if tags.has('build_pdflatex'):
   extensions.append('mcrl2_pdflatex')

if tags.has('build_manual'):
   extensions.append('mcrl2_manual')

# -- App setup - executed before the build process starts --------------------
def setup(app):
    os.makedirs(_SPHINX_BUILD_TEMP_DIR, mode = 0o755, exist_ok = True)
    os.makedirs(_SPHINX_BUILD_OUT_DIR, mode = 0o755, exist_ok = True)

    # Copy mcrl2_package.sty to the build temp dir (a path without spaces) and
    # add it to TEXINPUTS. LaTeX's \usepackage cannot handle paths with spaces,
    # so we avoid embedding the source path in imgmath_latex_preamble.
    shutil.copy2(Path(__file__).parent / 'mcrl2_package.sty', _SPHINX_BUILD_TEMP_DIR)
    tex_inputs = os.environ.get('TEXINPUTS', '')
    os.environ['TEXINPUTS'] = f'{_SPHINX_BUILD_TEMP_DIR}:{tex_inputs}' if tex_inputs else f'{_SPHINX_BUILD_TEMP_DIR}::'

    olddir = os.getcwd()
    try:
        os.chdir(_SPHINX_BUILD_TEMP_DIR)
    finally:
        os.chdir(olddir)
