# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# import sphinx_rtd_theme
# sys.path.insert(0, os.path.abspath('.'))
# sys.path.insert(0, os.path.abspath('..'))
# sys.path.append(os.path.dirname(__file__))


# -- Project information -----------------------------------------------------

project = 'DICAST'
copyright = '2021'
author = 'A. Fenn, O. Tsoy, A. Dietrich, T. Faro, F. Rößler'

# The full version, including alpha/beta/rc tags
release = '1.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'recommonmark',
    'sphinx_rtd_theme',
    'sphinx-prompt'
#    'sphinx.ext.autosectionlabel'
]

#needed for readthedocs.in
master_doc = 'contents'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Git
gitlab_url = 'ge46ban/dockers'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'collapse_navigation' : False,
    'sticky_navigation' : True,
    'navigation_depth' : 3,
    'titles_only' : True,
     'logo_only': True

}
html_logo = 'img/dicast_logo.svg'
html_favicon = 'img/favicon.ico'
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
# custom.css is inside one of the html_static_path folders (e.g. _static)
html_css_files = ["_templates/custom.css"]

def setup(app):
    app.add_css_file('_templates/custom.css')

# def setup(app):
#     app.add_css_file('custom.css')

# sphinx-notfound-page
# https://github.com/readthedocs/sphinx-notfound-page
notfound_context = {
    'title': 'Page Not Found',
    'body': '''
<h1>Page Not Found</h1>
<p>Sorry, we couldn't find that page.</p>
<p>Try using the search box or go to the homepage.</p>
''',
}