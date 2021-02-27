"""
dynophores.core.dynophore

Handles dynophore class.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from dynophores.core.superfeature import SuperFeature
from dynophores.core.envpartner import EnvPartner


class Dynophore:
    """
    Class to store dynophore data, i.e. data on superfeatures and their environmental partners,
    and important functions to interact with this data.

    Attributes
    ----------
    id : str
        Dynophore name.
    superfeatures : list of dynophores.base.SuperFeature
        Dynophore superfeatures.
    """

    def __init__(self):

        self.id = None
        self.superfeatures = []

    @classmethod
    def from_files(cls, dynophore_path):
        """
        Load dynophore data from DynophoreApp directory as Dynophore instance.

        Parameters
        ----------
        dynophore_path : pathlib.Path
            Path to DynophoreApp folder.

        Returns
        -------
        dynophores.Dynophore
            Dynophore.
        """

        dynophore = cls()

        # Get all files and filename components
        dynophore_path_options = [Path(dynophore_path) / "data", Path(dynophore_path) / "raw_data"]
        if dynophore_path_options[0].exists():
            dynophore_files = [file for file in dynophore_path_options[0].glob("*")]
            dynophore_files_components = [dynophore._file_components(i) for i in dynophore_files]
            print(f"Read files from {dynophore_path_options[0]}.")
        elif dynophore_path_options[1].exists():
            dynophore_files = [file for file in dynophore_path_options[1].glob("*")]
            dynophore_files_components = [
                dynophore._file_components_alternative(i) for i in dynophore_files
            ]
            print(f"Read files from {dynophore_path_options[1]}.")
        else:
            raise FileNotFoundError(
                "Could not find directory: "
                f"{' or '.join([str(i) for i in dynophore_path_options])}"
            )

        # Iterate over superfeatures
        superfeatures = []
        for superfeature_file_components in [
            i for i in dynophore_files_components if i["envpartner_id"] is None
        ]:
            # Iterate over environmental partners
            envpartners = []
            for envpartner_file_components in [
                i
                for i in dynophore_files_components
                if i["superfeature_id"] == superfeature_file_components["superfeature_id"]
                and i["envpartner_id"] is not None
            ]:
                # Get distances and occurrences for environmental partner
                try:
                    distances = np.loadtxt(
                        fname=envpartner_file_components["filepath"],
                        dtype=int,
                        delimiter=",",
                        usecols=1,
                    )
                    occurrences = np.loadtxt(
                        fname=envpartner_file_components["filepath"],
                        dtype=float,
                        delimiter=",",
                        usecols=0,
                    )
                except IndexError:
                    # In case data comes from master thesis DynophoreApp output
                    distances = np.loadtxt(
                        fname=envpartner_file_components["filepath"],
                        dtype=int,
                        delimiter=" ",
                        usecols=1,
                    )
                    occurrences = np.loadtxt(
                        fname=envpartner_file_components["filepath"],
                        dtype=float,
                        delimiter=" ",
                        usecols=0,
                    )

                # Set environmental partner
                envpartner = EnvPartner(
                    envpartner_file_components["envpartner_id"],
                    envpartner_file_components["envpartner_residue_name"],
                    envpartner_file_components["envpartner_residue_number"],
                    envpartner_file_components["envpartner_chain"],
                    envpartner_file_components["envpartner_atom_numbers"],
                    distances,
                    occurrences,
                )
                envpartners.append(envpartner)

            # Sort environmental partners list by (first) atom number (alphabetically)
            envpartners = sorted(envpartners, key=lambda envpartner: envpartner.atom_numbers[0])

            # Set superfeature
            superfeature = SuperFeature(
                superfeature_file_components["superfeature_id"],
                superfeature_file_components["superfeature_feature_type"],
                superfeature_file_components["superfeature_atom_numbers"],
                np.loadtxt(fname=superfeature_file_components["filepath"], dtype=int),
                envpartners,
            )
            superfeatures.append(superfeature)

        # Sort superfeatures list by ID (alphabetically)
        superfeatures = sorted(superfeatures, key=lambda superfeature: superfeature.id)

        # Add superfeature to dynophore
        dynophore.id = dynophore_files_components[0]["dynophore_id"]
        dynophore.superfeatures = superfeatures

        return dynophore

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
            [superfeature.occurrences for superfeature in self.superfeatures],
            index=[superfeature.id for superfeature in self.superfeatures],
        ).transpose()

        # Sort columns by feature type (alphabetically)
        superfeature_ids = occurrence_superfeatures.columns.to_list()
        superfeature_ids.sort()
        occurrence_superfeatures = occurrence_superfeatures[superfeature_ids]

        return occurrence_superfeatures

    @property
    def envpartners_occurrences(self):
        """
        For each superfeature, get its environmental partners' occurrences per environmental
        partner and frame.

        Returns
        -------
        dict of pandas.DataFrame
            For each superfeature, occurrences (0=no, 1=yes) of an environmental partner (columns)
            in each frame (row).
        """

        return self._envpartners_data(type="occurrences")

    @property
    def envpartners_distances(self):
        """
        For each superfeature, get its environmental partners' distances per environmental partner
        and frame.

        Returns
        -------
        dict of pandas.DataFrame
            For each superfeature, distances to an environmental partner (columns) in each frame
            (row).
        """

        return self._envpartners_data(type="distances")

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

        return len(self.superfeatures[0].occurrences)

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
            {superfeature.id: superfeature.count for superfeature in self.superfeatures}
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

    @staticmethod
    def _file_components(filepath):
        """
        Get components from dynophore filename.

        Parameters
        ----------
        filepath : str or pathlib.Path
            Path to dynophore file.

        Returns
        -------
        dict
            Components for dynophore filename.
        """

        filepath = Path(filepath)

        file_split = filepath.stem.split("_")

        file_components = {
            "filepath": None,
            "dynophore_id": None,
            "superfeature_id": None,
            "superfeature_feature_type": None,
            "superfeature_atom_numbers": None,
            "envpartner_id": None,
            "envpartner_residue_name": None,
            "envpartner_residue_number": None,
            "envpartner_chain": None,
            "envpartner_atom_numbers": None,
        }

        # Example filepath
        # 1KE7-1_data_superfeature_H[4599,4602,4601,4608,4609,4600]_100.0.txt
        # Is split into
        # ['1KE7-1', 'data', 'superfeature', 'H[4599,4602,4601,4608,4609,4600]', '100.0']

        file_components["filepath"] = filepath
        file_components["dynophore_id"] = file_split[0]
        file_components["superfeature_id"] = file_split[3]
        file_components["superfeature_feature_type"] = file_components["superfeature_id"].split(
            "["
        )[0]
        file_components["superfeature_atom_numbers"] = [
            int(atom) for atom in file_components["superfeature_id"].split("[")[1][:-1].split(",")
        ]

        if len(file_split) == 10:

            # Example filepath
            # 1KE7-1_data_superfeature_HBA[4619]_12.3_envpartner_ASP_86_A[1313]_1.6.txt
            # Is split into
            # ['1KE7-1', 'data', 'superfeature', 'HBA[4619]', '12.3', 'envpartner', 'ASP', '86',
            # 'A[1313]', '1.6']

            file_components["envpartner_id"] = "-".join(file_split[6:9])
            file_components["envpartner_residue_name"] = file_split[6]
            file_components["envpartner_residue_number"] = int(file_split[7])
            file_components["envpartner_chain"] = file_split[8].split("[")[0]
            file_components["envpartner_atom_numbers"] = [
                int(atom) for atom in file_split[8].split("[")[1][:-1].split(",")
            ]

        return file_components

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

    @staticmethod
    def _file_components_alternative(filepath):
        """
        Get components from dynophore filename that follows the syntax used in the master thesis
        version of the DynophoreApp (from 2015).

        Parameters
        ----------
        filepath : str or pathlib.Path
            Path to dynophore file.

        Returns
        -------
        dict
            Components for dynophore filename.
        """

        filepath = Path(filepath)

        file_split = filepath.stem.split("_")

        file_components = {
            "filepath": None,
            "dynophore_id": None,
            "superfeature_id": None,
            "superfeature_feature_type": None,
            "superfeature_atom_numbers": None,
            "envpartner_id": None,
            "envpartner_residue_name": None,
            "envpartner_residue_number": None,
            "envpartner_chain": None,
            "envpartner_atom_numbers": None,
        }

        # Example filepath
        # 1KE7-1_data_superfeature_H[4599,4602,4601,4608,4609,4600]_100.0.txt
        # Is split into
        # ['1KE7-1', 'data', 'superfeature', 'H[4599,4602,4601,4608,4609,4600]', '100.0']

        file_components["filepath"] = filepath
        file_components["dynophore_id"] = file_split[0]
        file_components["superfeature_id"] = file_split[3].split("%")[0]
        file_components["superfeature_feature_type"] = file_components["superfeature_id"].split(
            "["
        )[0]
        file_components["superfeature_atom_numbers"] = [
            int(atom) for atom in file_components["superfeature_id"].split("[")[1][:-1].split(",")
        ]

        if len(file_split) == 8:

            # Example filepath
            # 1KE7-1_data_superfeature_HBA[4619]_12.3_envpartner_ASP_86_A[1313]_1.6.txt
            # Is split into
            # ['1KE7-1', 'data', 'superfeature', 'HBA[4619]', '12.3', 'envpartner', 'ASP', '86',
            # 'A[1313]', '1.6']

            file_components["envpartner_id"] = "-".join(file_split[5:8]).split("%")[0]
            file_components["envpartner_residue_name"] = file_split[5]
            file_components["envpartner_residue_number"] = int(file_split[6])
            file_components["envpartner_chain"] = file_split[7].split("%")[0].split("[")[0]
            file_components["envpartner_atom_numbers"] = [
                int(atom) for atom in file_split[7].split("%")[0].split("[")[1][:-1].split(",")
            ]

        return file_components

    def _envpartners_data(self, type="occurrences"):
        """
        Get occurrences or distances of all superfeatures' environmental partners.

        Parameters
        ----------
        type : str
            Data type: occurrences (default) or distances.

        Returns
        -------
        dict of DataFrame
            Occurrences (default) or distances for a superfeature's (dict key) environmental
            partners (columns) for all frames (rows).
        """

        types = ["occurrences", "distances"]
        if type in types:
            pass
        else:
            raise KeyError(f'Wrong type. Select from: {", ".join(types)}')

        envpartners = {}

        for superfeature in self.superfeatures:
            superfeature_envpartners = pd.DataFrame(
                [getattr(envpartner, type) for envpartner in superfeature.envpartners],
                index=[envpartner.id for envpartner in superfeature.envpartners],
            ).transpose()

            envpartners[superfeature.id] = superfeature_envpartners

        return envpartners
