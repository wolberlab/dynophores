"""
Dynophores
Dynamic pharmacophore modeling of molecular interactions
"""

# Add imports here
from .core.dynophore import Dynophore
from .viz import plot
from .viz import view3d

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
