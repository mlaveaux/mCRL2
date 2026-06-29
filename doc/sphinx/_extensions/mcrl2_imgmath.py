"""Fix for sphinx.ext.imgmath parallel-write race condition.

When sphinx builds with -j > 1, worker processes share a single imgmath tempdir
and race to write math.tex / math.dvi.  A partially-overwritten DVI file causes:

    DVI error: invalid DVI file (missing preamble) at position -1

Fix: replace the imgmath math-renderer registry entries with thin wrappers
defined in *this* module.  Sphinx pickles those function references as
``mcrl2_imgmath.safe_html_visit_math`` / ``safe_html_visit_displaymath``.
When a worker reconstructs them it imports this module, which at module level
patches ``sphinx.ext.imgmath.ensure_tempdir`` to return a per-process temporary
directory, eliminating the race.
"""

import os
import tempfile

import sphinx.ext.imgmath as _imgmath
from sphinx.ext.imgmath import (
    html_visit_math as _orig_visit_math,
    html_visit_displaymath as _orig_visit_displaymath,
)


# Runs in every process that imports this module (main + each worker).
def _per_process_tempdir(builder):
    attr = f'_imgmath_tempdir_{os.getpid()}'
    if not hasattr(builder, attr):
        setattr(builder, attr, tempfile.mkdtemp())
    return getattr(builder, attr)


_imgmath.ensure_tempdir = _per_process_tempdir


def safe_html_visit_math(self, node):
    _orig_visit_math(self, node)


def safe_html_visit_displaymath(self, node):
    _orig_visit_displaymath(self, node)


def _patch_registry(app):
    reg = app.registry
    if 'imgmath' in reg.html_inline_math_renderers:
        reg.html_inline_math_renderers['imgmath'] = (safe_html_visit_math, None)
    if 'imgmath' in reg.html_block_math_renderers:
        reg.html_block_math_renderers['imgmath'] = (safe_html_visit_displaymath, None)


def setup(app):
    app.connect('builder-inited', _patch_registry)
    return {'parallel_read_safe': True, 'parallel_write_safe': True}
