"""
Dynophores
Dynamic pharmacophore modeling of molecular interactions
"""

# Add imports here
from .core.dynophore import Dynophore
from .core.superfeature import SuperFeature
from .core.envpartner import EnvPartner
from .viz import plot

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
