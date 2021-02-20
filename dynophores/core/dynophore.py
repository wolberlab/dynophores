"""
dynophores.core.dynophore

Handles dynophore class.
"""

import numpy as np
import pandas as pd

from dynophores.core.superfeature import Superfeature
from dynophores.core.envpartner import EnvPartner


class Dynophore:
    """
    Class to store dynophore data, i.e. data on superfeatures and their environmental partners,
    and important functions to interact with this data.

    Attributes
    ----------
    id : str
        Dynophore name.
    superfeatures : list of dynophores.base.Superfeature
        Dynophore superfeatures.
    """

    def __init__(self):

        self.id = None
        self.superfeatures = []

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

        return self._get_envpartners_data(type="occurrences")

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

        return self._get_envpartners_data(type="distances")

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

        dynophore_count = pd.DataFrame(
            {superfeature.id: superfeature.count for superfeature in self.superfeatures}
        )

        return dynophore_count.apply(lambda x: round(x / self.n_frames * 100, 2))

    def is_superfeature(self, superfeature_name):
        """
        Check if a superfeature is part of the dynophore (by name).

        Parameters
        ----------
        superfeature_name : str
            Superfeature name
        """

        if superfeature_name not in self.envpartners_occurrences.keys():
            raise KeyError(f"Superfeature name {superfeature_name} is unknown.")

    def from_file(self, dynophore_path):
        """
        Load dynophore data from file to Dynophore instance.

        Parameters
        ----------
        dynophore_path : pathlib.Path
            Path to folder with dynophore data.
        """

        # Get all files and filename components
        dynophore_files = [file for file in dynophore_path.glob("*")]
        dynophore_files_components = [self._get_file_components(i) for i in dynophore_files]

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
                # Set environmental partner
                envpartner = EnvPartner(
                    envpartner_file_components["envpartner_id"],
                    envpartner_file_components["envpartner_residue_name"],
                    envpartner_file_components["envpartner_residue_number"],
                    envpartner_file_components["envpartner_chain"],
                    envpartner_file_components["envpartner_atom_numbers"],
                    np.loadtxt(
                        fname=envpartner_file_components["filepath"],
                        dtype=int,
                        delimiter=",",
                        usecols=1,
                    ),
                    np.loadtxt(
                        fname=envpartner_file_components["filepath"],
                        dtype=float,
                        delimiter=",",
                        usecols=0,
                    ),
                )
                envpartners.append(envpartner)

            # Set superfeature
            superfeature = Superfeature(
                superfeature_file_components["superfeature_id"],
                superfeature_file_components["superfeature_feature_type"],
                superfeature_file_components["superfeature_atom_numbers"],
                np.loadtxt(fname=superfeature_file_components["filepath"], dtype=int),
                envpartners,
            )
            superfeatures.append(superfeature)

        # Add superfeature to dynophore
        self.id = dynophore_files_components[0]["dynophore_id"]
        self.superfeatures = superfeatures

    def show_2d_dynophore(self):
        """
        Show 2D dynophore representation.

        Returns
        -------
        TBA
        """

        raise NotImplementedError("Not yet implemented.")

    def show_3d_dynophore(self):
        """
        Show 3D dynophore representation (using nglview).

        Returns
        -------
        TBA
        """

        raise NotImplementedError("Not yet implemented.")

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

    def _get_envpartners_data(self, type="occurrences"):
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
            raise ValueError(f'Wrong type. Select from: {", ".join(types)}')

        envpartners = {}

        for superfeature in self.superfeatures:
            superfeature_envpartners = pd.DataFrame(
                [getattr(envpartner, type) for envpartner in superfeature.envpartners],
                index=[envpartner.id for envpartner in superfeature.envpartners],
            ).transpose()

            envpartners[superfeature.id] = superfeature_envpartners

        return envpartners
