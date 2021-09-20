r"""
Welcome to dynophores! Below you will find the submodules exposed by the API

.. autosummary::
   :toctree: autosummary
   :nosignatures:

    Dynophore
    plot
    view2d
    view3d
"""

# Add imports here
from .core.dynophore import Dynophore
from .viz import plot, view2d, view3d

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
