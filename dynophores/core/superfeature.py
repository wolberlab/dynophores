"""
Handles the SuperFeature class, which describes one superfeature for one dynophore.
"""

import pandas as pd

from dynophores.core.envpartner import EnvPartner
from dynophores.core.chemicalfeaturecloud3d import ChemicalFeatureCloud3D


class SuperFeature:
    """
    Class to store superfeature data, i.e. data on environmental partners, and important functions
    to interact with this data.

    Attributes
    ----------
    id : str
        Superfeature ID.
    envpartner_ids : list of str
        Environmental partner identifiers available for this superfeature.
    feature_type : str
        Pharmacophoric feature type.
    atom_numbers : list of int
        List of atom IDs.
    occurrences : np.array
        Occurrence of superfeature (0=no, 1=yes) in each frame.
    envpartners : list of EnvPartner
        Superfeature's environmental partners.
    color : str
        Superfeature color.
    cloud : ChemicalFeatureCloud3D
        Chemical feature cloud.
    """

    def __init__(self, id, feature_type, atom_numbers, occurrences, envpartners, color, cloud):
        self.id = id
        self.feature_type = feature_type
        self.atom_numbers = atom_numbers
        self.occurrences = occurrences
        self.envpartners = {
            envpartner_id: envpartner
            if isinstance(envpartner, EnvPartner)
            else EnvPartner(**envpartner)
            for envpartner_id, envpartner in envpartners.items()
        }
        self.color = color
        self.cloud = (
            cloud if isinstance(cloud, ChemicalFeatureCloud3D) else ChemicalFeatureCloud3D(**cloud)
        )

    @property
    def envpartners_occurrences(self):
        """
        Get the superfeature's environmental partners' occurrences per environmental partner and
        frame.

        Returns
        -------
        pandas.DataFrame
            Occurrences (0=no, 1=yes) of an environmental partner (columns) in each frame (row).
        """

        return self._envpartners_occurrences(self._data(type="occurrences"))

    @property
    def envpartners_occurrences_collapsed(self):
        """
        Get the superfeature's environmental partners' occurrences per environmental partner and
        frame.
        If an environmental partner interacts multiple times with the same superfeature,
        aggregate them. This can happen if differen atoms of an environmental partner are involved
        in the same superfeature.

        Returns
        -------
        pandas.DataFrame
            Occurrences (0=no, 1=yes) of an environmental partner (columns) in each frame (row).
        """

        return self._envpartners_occurrences(self._data_collapsed())

    def _envpartners_occurrences(self, method_data):
        """
        Get the superfeature's environmental partners' occurrences per environmental partner and
        frame.

        Returns
        -------
        pandas.DataFrame
            Occurrences (0=no, 1=yes) of an environmental partner (columns) in each frame (row).
        """

        occurrences = method_data.astype("int32")

        # Sort columns by superfeature occurrence
        sorted_columns = occurrences.sum().sort_values(ascending=False).index
        occurrences = occurrences[sorted_columns]

        return occurrences

    @property
    def envpartners_distances(self):
        """
        Get the superfeature's environmental partners' distances per environmental partner and
        frame.

        Returns
        -------
        pandas.DataFrame
            Distances to an environmental partner (columns) in each frame (row)
        """

        distances = self._data(type="distances")

        # Sort columns by superfeature occurrence
        distances = distances[self.envpartners_occurrences.columns]

        return distances

    @property
    def n_frames(self):
        """
        Get superfeatures's number of frames.

        Returns
        -------
        int
            Number of frames.
        """

        return len(self.occurrences)

    @property
    def count(self):
        """
        Get number of frames in which the superfeature occurs, including the superfeature's
        environmental partners occurrences.

        Returns
        -------
        pandas.Series
            Superfeature count: The Series shows interactions (yes/no) to each single
            environmental partner as well as any environmental partner.
        """

        return self._count(self.envpartners_occurrences)

    @property
    def count_collapsed(self):
        """
        Get number of frames in which the superfeature occurs, including the superfeature's
        environmental partners occurrences (collapsed if they share the same residue!).

        Returns
        -------
        pandas.Series
            Superfeature count: The Series shows interactions (yes/no) to each single
            environmental partner as well as any environmental partner.
        """

        return self._count(self.envpartners_occurrences_collapsed)

    def _count(self, property_envpartners_occurrences):
        """
        Count the occurrence of the superfeature's environmental partners.

        Parameter
        ---------
        property : property_envpartners_occurrences
            If you want un-collapsed environmental partners, use `self.envpartners_occurrences`.
            If you want collapsed environmental partners, use
            `self.envpartners_occurrences_collapsed`.

        Returns
        -------
        pandas.Series
            Superfeature count: The Series shows interactions (yes/no) to each single
            environmental partner as well as any environmental partner.
        """

        superfeature_count = pd.Series(
            {"any": (property_envpartners_occurrences.sum(axis=1) != 0).sum()}
        )
        envpartners_count = property_envpartners_occurrences.sum()

        return superfeature_count.append(envpartners_count)

    @property
    def frequency(self):
        """
        Get frequency of frames in which the superfeature occurs, including the superfeature's
        environmental partners occurrences.

        Returns
        -------
        pandas.Series
            Superfeature frequency: The Series shows interactions (yes/no) to each single
            environmental partner as well as any environmental partner.
        """

        return self._frequency(self.count)

    @property
    def frequency_collapsed(self):
        """
        Get frequency of frames in which the superfeature occurs, including the superfeature's
        environmental partners occurrences (collapsed if they share the same residue!).

        Returns
        -------
        pandas.Series
            Superfeature frequency: The Series shows interactions (yes/no) to each single
            environmental partner as well as any environmental partner.
        """

        return self._frequency(self.count_collapsed)

    def _frequency(self, property_count):
        """
        Get the frequency of the occurrence of the superfeature's environmental partners.

        Parameter
        ---------
        property : property_count
            If you want un-collapsed environmental partners, use `self.count`.
            If you want collapsed environmental partners, use `self.count_collapsed`.

        Returns
        -------
        pandas.Series
            Superfeature frequency: The Series shows interactions (yes/no) to each single
            environmental partner as well as any environmental partner.
        """

        return property_count.apply(lambda x: round(x / self.n_frames * 100, 2))

    def _data(self, type="occurrences"):
        """
        Get occurrences or distances of a superfeature's environmental partners.

        Parameters
        ----------
        type : str
            Data type: occurrences (default) or distances.

        Returns
        -------
        DataFrame
            Occurrences (default) or distances of a superfeature's environmental partners
            (columns) for all frames (rows).
        """

        types = ["occurrences", "distances"]
        if type in types:
            pass
        else:
            raise KeyError(f'Wrong type. Select from: {", ".join(types)}')

        return pd.DataFrame(
            {
                envpartner_id: getattr(envpartner, type)
                for envpartner_id, envpartner in self.envpartners.items()
            }
        )

    def _data_collapsed(self):
        """TODO"""

        # List of environmental partner IDs (e.g. ILE-10-A[169,171,172])
        ids = self.envpartners_occurrences.columns
        # Unique list of residue IDs (e.g. ILE-10-A)
        residue_ids = [envpartner.residue_id for _, envpartner in self.envpartners.items()]
        residue_ids = list(set(residue_ids))

        occurrences_dict = {}

        # For each unique residue ID,
        # we want to aggregate data for all environmental partners that belong to the same residue
        for residue_id in residue_ids:
            # Get all environmental partner IDs that belong to this residue
            ids_to_be_collapsed = [_id for _id in ids if _id.startswith(residue_id)]

            # Merge all atom numbers
            atom_numbers = []
            for _id in ids_to_be_collapsed:
                atom_numbers.extend(self.envpartners[_id].atom_numbers)
            atom_numbers = sorted(list(set(atom_numbers)))
            id_collapsed = f"{residue_id}{atom_numbers}".replace(" ", "")

            # Merge all occurrences
            occurrences = [self.envpartners[_id].occurrences for _id in ids_to_be_collapsed]
            occurrences = pd.DataFrame(occurrences)
            # If frame 1 in any environmental partner, set to 1 in collapsed environmental partner
            occurrences = occurrences.sum().apply(lambda x: 1 if x > 0 else 0)

            occurrences_dict[id_collapsed] = occurrences

        occurrences = pd.DataFrame(occurrences_dict, dtype="int32")

        return occurrences
