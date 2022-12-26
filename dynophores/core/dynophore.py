"""
Handles the Dynophore class, describing a dynophore with its superfeatures and environmental
partners.
"""

from pathlib import Path
import logging

import pandas as pd

from dynophores import parsers
from dynophores.core.superfeature import SuperFeature
from dynophores.core.ligand import Ligand

logger = logging.getLogger(__name__)


class Dynophore:
    """
    Class to store dynophore data, i.e. data on superfeatures and their environmental partners,
    and important functions to interact with this data.

    Attributes
    ----------
    id : str
        Dynophore name.
    ligand : dynophores.Ligand
        Ligand information.
    superfeatures : dict of str: dynophores.SuperFeature
        Dynophore superfeatures (values) by superfeature IDs (keys).
    """

    def __init__(
        self,
        id,
        ligand,
        superfeatures,
        **kwargs,
    ):

        self.id = id
        self.ligand = ligand if isinstance(ligand, Ligand) else Ligand(**ligand)
        self.superfeatures = {
            superfeature_id: superfeature
            if isinstance(superfeature, SuperFeature)
            else SuperFeature(**superfeature)
            for superfeature_id, superfeature in superfeatures.items()
        }

    @classmethod
    def from_dir(cls, dynophore_path):
        """
        Load dynophore data from DynophoreApp output directory.

        Parameters
        ----------
        dynophore_path : pathlib.Path
            Path to DynophoreApp output folder.

        Returns
        -------
        dynophores.Dynophore
            Dynophore.
        """

        dynophore_path = Path(dynophore_path)

        if dynophore_path.is_dir():

            # Set JSON path
            json_path = list(dynophore_path.glob("*.json"))
            if len(json_path) == 1:
                json_path = json_path[0]
            else:
                raise ValueError(
                    f"None or too many JSON files in {dynophore_path}. Only one allowed."
                )

            # Set PML path
            pml_path = list(dynophore_path.glob("*.pml"))
            if len(pml_path) == 1:
                pml_path = pml_path[0]
            else:
                raise ValueError(
                    f"None or roo many PML files in {dynophore_path}. Only one allowed."
                )

            # Create Dynophore object
            return cls.from_files(json_path, pml_path)

        else:
            raise FileNotFoundError(
                "Input directory does not exist or is no directory. "
                f"Your input: {dynophore_path}"
            )

    @classmethod
    def from_files(cls, json_path, pml_path):
        """
        Load dynophore data from JSON and PML file.

        Parameters
        ----------
        json_path : pathlib.Path
            Path to dynophore JSON file.
        pml_path : pathlib.Path
            Path to dynophore PML file.

        Returns
        -------
        dynophores.Dynophore
            Dynophore.
        """

        dynophore_dict = parsers._json_pml_to_dict(json_path, pml_path)
        dynophore = cls(**dynophore_dict)
        return dynophore

    @property
    def clouds(self):
        """
        Dynophore clouds.

        Returns
        -------
        dict of pandas.DataFrame
            Per superfeature, cloud point coordinates.
        """

        return {
            superfeature_id: superfeature.cloud.data
            for superfeature_id, superfeature in self.superfeatures.items()
        }

    def cloud_by_superfeature(self, superfeature_id):
        """
        Superfeature cloud.

        Parameters
        ----------
        superfeature_id : str
            Superfeature identifier.

        Returns
        -------
        pandas.DataFrame
            Superfeature cloud point coordinates.
        """

        return self.superfeatures[superfeature_id].cloud.data

    @property
    def superfeatures_occurrences(self):
        """
        Get the superfeatures' occurrences per superfeature and frame.

        Returns
        -------
        pandas.DataFrame
            Occurrences (0=no, 1=yes) of a superfeature (columns) in each frame (row).
        """

        occurrences = pd.DataFrame(
            {
                superfeature_id: superfeature.occurrences
                for superfeature_id, superfeature in self.superfeatures.items()
            }
        ).astype("int32")

        # Sort columns by superfeature occurrence
        sorted_columns = occurrences.sum().sort_values(ascending=False).index
        occurrences = occurrences[sorted_columns]

        return occurrences

    def envpartners_occurrences_by_superfeature(self, superfeature_id):
        """
        For a given superfeature, get its environmental partners' occurrences per environmental
        partner and frame.

        Parameters
        ----------
        superfeature_id : str
            Superfeature identifier.

        Returns
        -------
        pandas.DataFrame
            For a given superfeature, occurrences (0=no, 1=yes) of an environmental partner
            (columns) in each frame (row).
        """

        superfeature = self.superfeatures[superfeature_id]
        return superfeature.envpartners_occurrences

    @property
    def envpartners_occurrences(self):
        """
        For each superfeature, get its environmental partners' occurrences per environmental
        partner and frame.

        Returns
        -------
        dict of pandas.DataFrame
            For each superfeature (keys), occurrences (0=no, 1=yes) of an environmental partner
            (columns) in each frame (row).
        """

        return {
            superfeature_id: superfeature.envpartners_occurrences
            for superfeature_id, superfeature in self.superfeatures.items()
        }

    def envpartners_distances_by_superfeature(self, superfeature_id):
        """
        For a given superfeature, get its environmental partners' distances per environmental
        partner and frame.

        Parameters
        ----------
        superfeature_id : str
            Superfeature identifier.

        Returns
        -------
        pandas.DataFrame
            For a given superfeature, distances to an environmental partner (columns) in each frame
            (row).
        """

        superfeature = self.superfeatures[superfeature_id]
        return superfeature.envpartners_distances

    @property
    def envpartners_distances(self):
        """
        For each superfeature, get its environmental partners' distances per environmental
        partner and frame.

        Returns
        -------
        dict of pandas.DataFrame
            For each superfeature (keys), distances to an environmental partner (columns) in each
            frame (row).
        """

        return {
            superfeature_id: superfeature.envpartners_distances
            for superfeature_id, superfeature in self.superfeatures.items()
        }

    @property
    def n_superfeatures(self):
        """
        Get dynophore's number of superfeatures.

        Returns
        -------
        int
            Number of superfeatures.
        """

        return len(self.superfeatures)

    @property
    def n_frames(self):
        """
        Get dynophore's number of frames.

        Returns
        -------
        int
            Number of frames.
        """

        superfeature = next(iter(self.superfeatures.values()))
        return len(superfeature.occurrences)

    @property
    def count(self):
        """
        Get number of frames in which each dynophore occurs, including the superfeatures and
        superfeatures' environmental partners occurrences.

        Returns
        -------
        pandas.DataFrame
            Dynophore count: The DataFrame shows interaction (yes/no) for superfeatures (rows) to
            each single environmental partner as well as any environmental partner (columns).
        """

        dynophore_count = pd.DataFrame(
            {
                superfeature_id: superfeature.count
                for superfeature_id, superfeature in self.superfeatures.items()
            }
        )
        dynophore_count.fillna(0, inplace=True)
        dynophore_count = dynophore_count.astype("int32")

        return dynophore_count

    @property
    def frequency(self):
        """
        Get frequency of frames in which each dynophore occurs, including the superfeatures and
        superfeatures' environmental partners occurrences.

        Returns
        -------
        pandas.DataFrame
            Dynophore frequency: The DataFrame shows interaction (yes/no) for superfeatures (rows)
            to each single environmental partner as well as any environmental partner (columns).
        """

        return self.count.apply(lambda x: round(x / self.n_frames * 100, 2))

    @property
    def unique_envpartners_chain_residue_number(self):
        """
        List of unique environmental partners (chain and residue number).
        Useful for 3D visualization of interacting pocket residues.

        Returns
        -------
        list of tuple (str, int)
            List of (chain, residue number) tuples.
        """

        envpartners = [
            (envpartner.chain, envpartner.residue_number)
            for superfeature_id, superfeature in self.superfeatures.items()
            for envpartner_id, envpartner in superfeature.envpartners.items()
        ]
        envpartners = set(envpartners)
        return envpartners

    def _raise_keyerror_if_invalid_superfeature_id(self, superfeature_id):
        """
        Check if dynophore has a certain superfeature (by name).

        Parameters
        ----------
        superfeature_id : str
            Superfeature ID
        """

        if superfeature_id not in self.envpartners_occurrences.keys():
            raise KeyError(f"Superfeature ID {superfeature_id} is unknown.")

    def _superfeature_ids_frequencies_strings(self, superfeature_ids):
        """
        Get superfeature IDs with frequencies as strings (useful for plotting).

        Parameters
        ----------
        superfeature_ids : list of str
            Superfeature IDs.

        Returns
        -------
        list of str
            Superfeature IDs with frequencies.
        """

        superfeature_ids_frequencies = round(self.frequency.loc["any", superfeature_ids], 1)
        superfeature_ids_frequencies_strings = (
            superfeature_ids_frequencies.reset_index()
            .apply(lambda x: f"{x['index']} {x['any']}%", axis=1)
            .to_list()
        )
        return superfeature_ids_frequencies_strings

    def _envpartner_names_frequencies_strings(self, superfeature_id, envpartners_names):
        """
        Get enpartner names for a superfeature with frequencies as strings (useful for plotting).

        Parameters
        ----------
        superfeature_id : str
            Superfeature.
        envpartners_names : list of str
            List of envpartners involved in that superfeature.

        Returns
        -------
        list of str
            Envpartners with frequencies.
        """

        envpartners_occurrences = self.envpartners_occurrences[superfeature_id]

        envpartner_names_frequencies = round(
            envpartners_occurrences.sum() / envpartners_occurrences.shape[0] * 100, 1
        )
        envpartner_names_frequencies = envpartner_names_frequencies[envpartners_names]
        envpartner_names_frequencies_strings = (
            envpartner_names_frequencies.reset_index()
            .apply(lambda x: f"{x['index']} {x[0]}%", axis=1)
            .to_list()
        )

        return envpartner_names_frequencies_strings

    @property
    def superfeatures_atom_serials(self):
        """TODO"""

        return {
            superfeature_id: superfeature.atom_numbers
            for superfeature_id, superfeature in self.superfeatures.items()
        }

    @property
    def superfeatures_colors(self):
        """TODO"""

        return {
            superfeature_id: superfeature.color
            for superfeature_id, superfeature in self.superfeatures.items()
        }
