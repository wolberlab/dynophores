"""
Contains base functions/classes.
"""

from pathlib import Path

import numpy as np
import pandas as pd


class Dynophore:
    """
    Class to store dynophore data, i.e. data on superfeatures and their environmental partners, and important functions
    to interact with this data.

    Attributes
    ----------
    id : str
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

        self.id = None
        self.superfeatures = []

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

        # Get all files
        dynophore_files = [file for file in dynophore_path.glob('*')]

        # Get filename components
        dynophore_files_components = [self._get_file_components(i) for i in dynophore_files]

        # Initialize dynophore
        self.id = dynophore_files_components[0]['dynophore_id']

        # Iterate over superfeatures
        for superfeature_file_components in [i for i in dynophore_files_components if i['envpartner_id'] is None]:

            # Initialize superfeature and set attributes
            superfeature = Superfeature()
            superfeature.id = superfeature_file_components['superfeature_id']
            superfeature.feature_type = superfeature_file_components['superfeature_feature_type']
            superfeature.atom_numbers = superfeature_file_components['superfeature_atom_numbers']
            superfeature.occurrences = np.loadtxt(fname=superfeature_file_components['filepath'], dtype=int)

            # Iterate over environmental partners
            for envpartner_file_components in [i for i in dynophore_files_components
                                               if i['superfeature_id'] == superfeature.id and
                                                  i['envpartner_id'] is not None]:
                # Initialize environmental partner and set attributes
                envpartner = EnvPartner()
                envpartner.id = envpartner_file_components['envpartner_id']
                envpartner.residue_name = envpartner_file_components['envpartner_residue_name']
                envpartner.residue_number = envpartner_file_components['envpartner_residue_number']
                envpartner.chain = envpartner_file_components['envpartner_chain']
                envpartner.atom_numbers = envpartner_file_components['envpartner_atom_numbers']
                envpartner.occurrences = np.loadtxt(fname=envpartner_file_components['filepath'], dtype=int,
                                                    delimiter=',', usecols=1)
                envpartner.distances = np.loadtxt(fname=envpartner_file_components['filepath'], dtype=float,
                                                  delimiter=',', usecols=0)

                # Add environmental partner to superfeature
                superfeature.envpartners.append(envpartner)

            # Add superfeature to dynophore
            self.superfeatures.append(superfeature)

    @staticmethod
    def _get_file_components(filepath):
        """
        Get components from dynophore filename.

        Parameters
        ----------
        filepath : pathlib.Path
            Path to dynophore file.

        Returns
        -------
        dict
            Components for dynophore filename.
        """

        file_split = filepath.stem.split('_')
        file_split.remove('')

        file_components = {
            'filepath': None,
            'dynophore_id': None,
            'superfeature_id': None,
            'superfeature_feature_type': None,
            'superfeature_atom_numbers': None,
            'envpartner_id': None,
            'envpartner_residue_name': None,
            'envpartner_residue_number': None,
            'envpartner_chain': None,
            'envpartner_atom_numbers': None
        }

        file_components['filepath'] = filepath
        file_components['dynophore_id'] = file_split[0]
        file_components['superfeature_id'] = file_split[3].split('%')[0]
        file_components['superfeature_feature_type'] = file_components['superfeature_id'].split('[')[0]
        file_components['superfeature_atom_numbers'] = [int(atom) for atom in
                                                        file_components['superfeature_id'].split('[')[1][:-1].split(
                                                            ',')]

        if len(file_split) == 6:
            file_components['envpartner_id'] = file_split[5].split('%')[0]
            file_components['envpartner_residue_name'] = file_components['envpartner_id'].split('-')[0]
            file_components['envpartner_residue_number'] = int(file_components['envpartner_id'].split('-')[1])
            file_components['envpartner_chain'] = file_components['envpartner_id'].split('-')[2].split('[')[0]
            file_components['envpartner_atom_numbers'] = [int(atom) for atom in
                                                          file_components['envpartner_id'].split('[')[1][:-1].split(
                                                              ',')]

        return file_components

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
    id : str
        Superfeature name.
    feature_type : str
        Pharmacophoric feature type.
    atom_numbers : list of int
        List of atom IDs.
    occurrences : pandas.DataFrame
        Per superfeature (columns), occurrences (0=no, 1=yes) in each frame (row).
    envpartners : list of EnvPartner
        Superfeature's environmental partners.
    """

    def __init__(self):

        self.id = None
        self.feature_type = None
        self.atom_numbers = None
        self.occurrences = None
        self.envpartners = []

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


class EnvPartner:
    """
    Class to store environmental partner data.

    Attributes
    ----------
    id : str
        ID of environmental partner.
    residue_name : str
        Residue name of Environmental partner.
    residue_number : int
        Residue number of environmental partner.
    chain : str
        Chain of environmental partner.
    atom_ids : list of int
        List of atom IDs.
    occurrences : numpy.array
        Occurrences of interaction with environmental partner (0=no, 1=yes) in each frame.
    distances : numpy.array
        Interaction distances in each frame.
    """

    def __init__(self):

        self.id = None
        self.residue_name = None
        self.residue_number = None
        self.chain = None
        self.atom_ids = None
        self.occurrences = None
        self.distances = None
