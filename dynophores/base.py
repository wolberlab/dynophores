"""
Contains base functions/classes.
"""

from pathlib import Path

import pandas as pd


class Dynophore:
    """
    Class to store dynophore data, i.e. data on superfeatures and their environmental partners, and important functions
    to interact with this data.

    Attributes
    ----------
    name : str
        Dynophore name.
    superfeatures : list of dynophores.base.Superfeature
        Dynophore superfeatures.

    Examples
    --------
    >>> from pathlib import Path
    >>> from dynophores.base import Dynophore
    >>> dynophore_path = Path(__name__).parent / 'dynophores' / 'tests' / 'data' / ''
    >>> dynophore = Dynophore()
    >>> dynophore.from_file(dynophore_path)
    """

    def __init__(self):

        self.name = None
        self.superfeatures = None

    @property
    def n_frames(self):
        """
        Get dynophore's number of frames.

        Returns
        -------
        int
            Number of frames.
        """

        pass

    @property
    def superfeatures_occurrences(self):
        """
        Get the superfeatures' occurrences per superfeature and frame.

        Returns
        -------
        pandas.DataFrame
            Occurrences (0=no, 1=yes) of a superfeature (columns) in each frame (row).
        """

        pass

    @property
    def envpartners_occurrences(self):
        """
        For each superfeature, get its environmental partners' occurrences per environmental partner and frame.

        Returns
        -------
        dict of pandas.DataFrame
            For each superfeature, occurrences (0=no, 1=yes) of an environmental partner (columns) in each frame (row).
        """

        pass

    @property
    def envpartners_distances(self):
        """
        For each superfeature, get its environmental partners' distances per environmental partner and frame.

        Returns
        -------
        dict of pandas.DataFrame
            For each superfeature, distances to an environmental partner (columns) in each frame (row).

        """

        pass

    @property
    def count(self, relative=False):
        """
        Get number of frames in which each superfeature occurs, including the superfeatures' environmental partners
        occurrences.

        Parameters
        ----------
        relative : bool
            Absolute (default) or relative count.

        Returns
        -------
        pandas.DataFrame
            Superfeatures count: The DataFrame shows interaction (yes/no) for superfeatures (rows) to each single
            environmental partner as well as any environmental partner (columns).
        """

        pass

    def from_file(self, dynophore_path):
        """
        Load dynophore data from file to Dynophore instance.

        Parameters
        ----------
        dynophore_path : pathlib.Path
            Path to folder with dynophore data.
        """

        pass

    def show_2d_dynophore(self):
        """
        Show 2D dynophore representation.

        Returns
        -------
        TBA
        """

        pass

    def show_3d_dynophore(self):
        """
        Show 3D dynophore representation (using nglview).

        Returns
        -------
        TBA
        """

        pass

    def plot_superfeatures_occurrences(self):
        """
        Plot the superfeatures' occurrences as barcode.

        Returns
        -------
        TBA
        """
        pass


class Superfeature:
    """
    Class to store superfeature data, i.e. data on environmental partners, and important functions to interact with this
    data.

    Attributes
    ----------
    name : str
        Superfeature name.
    feature_type : str
        Pharmacophoric feature type.
    atom_numbers : list of int
        List of atom numbers.
    occurrences : pandas.DataFrame
        Per superfeature (columns), occurrences (0=no, 1=yes) in each frame (row).
    envpartners : list of EnvPartner
        Superfeature's environmental partners.
    """

    def __init__(self):

        self.name = None
        self.feature_type = None
        self.atom_numbers = None
        self.occurrences = None
        self.envpartners = None

    @property
    def envpartners_occurrences(self):
        """
        Get the superfeature's environmental partners' occurrences per environmental partner and frame.

        Returns
        -------
        pandas.DataFrame
            Occurrences (0=no, 1=yes) of an environmental partner (columns) in each frame (row).

        """

        pass

    @property
    def envpartners_distances(self):
        """
        Get the superfeature's environmental partners' distances per environmental partner and frame.

        Returns
        -------
        pandas.DataFrame
            Distances to an environmental partner (columns) in each frame (row)
        """

        pass

    @property
    def count(self, relative=False):
        """
        Get number of frames in which the superfeature occurs, including the superfeature's environmental partners
        occurrences.

        Parameters
        ----------
        relative : bool
            Absolute (default) or relative count.

        Returns
        -------
        pandas.Series
            Superfeature count: The Series shows interactions (yes/no) to each single environmental partner as well as
            any environmental partner.
        """

        pass

    def from_file(self, superfeature_path, envpartners_paths):
        """
        Load superfeature data from files to Superfeature instance.

        Parameters
        ----------
        superfeature_path : pathlib.Path
            Path to file with superfeature data.
        envpartners_paths : list of pathlib.Path
            List of paths to files with superfeature's environmental partners data.
        """

        pass

    def plot_envpartners_occurrences(self):
        """
        Plot the superfeature's environmental occurrences as barcode.

        Returns
        -------
        TBA
        """

        pass

    def plot_envpartners_distances(self):
        """
        Plot the superfeature's environmental distances as series or histogram.

        Returns
        -------
        TBA
        """

        pass

