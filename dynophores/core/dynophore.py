"""
dynophores.core.dynophore

Handles the Dynophore class, describing a dynophore with its superfeatures and environmental
partners.
"""

from pathlib import Path
import logging

import numpy as np
import pandas as pd

from dynophores import parsers
from dynophores.core.superfeature import SuperFeature
from dynophores.core.envpartner import EnvPartner
from dynophores.core.chemicalfeaturecloud3d import ChemicalFeatureCloud3D


logger = logging.getLogger(__name__)


class Dynophore:
    """
    Class to store dynophore data, i.e. data on superfeatures and their environmental partners,
    and important functions to interact with this data.

    Attributes
    ----------
    id : str
        Dynophore name.
    superfeature_ids : list of str
        Superfeature identifiers available for this dynophore.
    superfeatures : list of dynophores.base.SuperFeature
        Dynophore superfeatures.
    """

    def __init__(self):

        self.id = None
        self.superfeature_ids = []
        self.superfeatures = []

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

        json_path = Path(json_path)
        pml_path = Path(pml_path)

        json_dict = parsers._json_to_dict(json_path)
        pml_dict = parsers._pml_to_dict(pml_path)

        # Check if superfeatures are the same in both files
        json_superfeatures = sorted([sp["id"] for sp in json_dict["superfeatures"]])
        pml_superfeatures = sorted(list(pml_dict.keys()))
        if json_superfeatures != pml_superfeatures:
            raise ValueError(
                f"Your PML and JSON files are not matching. Superfeatures must be the same.\n"
                f"Your JSON file: {json_path}\n"
                f"Your PML file: {pml_path}"
            )

        dynophore = cls()

        json_dynophore_dict = json_dict

        dynophore.id = json_dynophore_dict["id"]

        superfeatures = {}
        for json_superfeature_dict in json_dynophore_dict["superfeatures"]:

            envpartners = {}
            for json_envpartner_dict in json_superfeature_dict["envpartners"]:
                try:
                    residue_name = json_envpartner_dict["residue_name"]
                    residue_number = json_envpartner_dict["residue_number"]
                    chain = json_envpartner_dict["chain"]
                except KeyError:
                    residue_name = json_envpartner_dict["name"].split("_")[0]
                    residue_number = json_envpartner_dict["name"].split("_")[1]
                    chain = json_envpartner_dict["name"].split("_")[2]
                envpartner_id = json_envpartner_dict["id"]
                # TODO Remove next line, once updated in DynophoreApp
                envpartner_id = envpartner_id.replace("_", "-")
                envpartner = EnvPartner(
                    envpartner_id,
                    residue_name,
                    residue_number,
                    chain,
                    json_envpartner_dict["atom_numbers"],
                    np.array(json_envpartner_dict["occurrences"]),
                    np.array(json_envpartner_dict["distances"]),
                )
                envpartners[envpartner_id] = envpartner

            superfeature_id = json_superfeature_dict["id"]
            cloud = ChemicalFeatureCloud3D(**pml_dict[superfeature_id])

            superfeature = SuperFeature(
                superfeature_id,
                json_superfeature_dict["feature_type"],
                json_superfeature_dict["atom_numbers"],
                np.array(json_superfeature_dict["occurrences"]),
                envpartners,
                cloud,
            )
            superfeatures[superfeature_id] = superfeature
        dynophore.superfeatures = superfeatures
        dynophore.superfeature_ids = list(superfeatures.keys())

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

        occurrence_superfeatures = pd.DataFrame(
            {
                superfeature_id: superfeature.occurrences
                for superfeature_id, superfeature in self.superfeatures.items()
            }
        )

        # Sort columns by feature type (alphabetically)
        superfeature_ids = occurrence_superfeatures.columns.to_list()
        superfeature_ids.sort()
        occurrence_superfeatures = occurrence_superfeatures[superfeature_ids]

        return occurrence_superfeatures

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
            For a given superfeature, occurrences (0=no, 1=yes) of an environmental partner (columns)
            in each frame (row).
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

    def _raise_keyerror_if_invalid_superfeature_name(self, superfeature_name):
        """
        Check if dynophore has a certain superfeature (by name).

        Parameters
        ----------
        superfeature_name : str
            Superfeature name
        """

        if superfeature_name not in self.envpartners_occurrences.keys():
            raise KeyError(f"Superfeature name {superfeature_name} is unknown.")

    def _superfeature_names_frequencies_strings(self, superfeature_names):
        """
        Get superfeature names with frequencies as strings (useful for plotting).
        TODO unit test

        Parameters
        ----------
        superfeature_names : list of str
            Superfeature names.

        Returns
        -------
        list of str
            Superfeature names with frequencies.
        """

        superfeature_names_frequencies = round(self.frequency.loc["any", superfeature_names], 1)
        superfeature_names_frequencies_strings = (
            superfeature_names_frequencies.reset_index()
            .apply(lambda x: f"{x['index']} {x['any']}%", axis=1)
            .to_list()
        )
        return superfeature_names_frequencies_strings

    def _envpartner_names_frequencies_strings(self, superfeature_name, envpartners_names):
        """
        TODO docstring + unit test
        """

        envpartners_occurrences = self.envpartners_occurrences[superfeature_name]

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
