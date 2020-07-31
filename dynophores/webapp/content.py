"""
This module defines the webapp content.
"""

import logging

from .utils import EXAMPLE_PATH
from dynophores.base import Dynophore

logger = logging.getLogger(__name__)


def load_dynophore():

    dynophore = Dynophore()
    dynophore.from_file(EXAMPLE_PATH)

    return dynophore
