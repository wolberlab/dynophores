"""
dynophores.core.superfeature

Handles superfeature class.
"""

import pandas as pd


class SuperFeature:
    """
    Class to store superfeature data, i.e. data on environmental partners, and important functions
    to interact with this data.

    Attributes
    ----------
    id : str
        Superfeature name.
    feature_type : str
        Pharmacophoric feature type.
    atom_numbers : list of int
        List of atom IDs.
    occurrences : np.array
        Occurrence of superfeature (0=no, 1=yes) in each frame.
    envpartners : list of EnvPartner
        Superfeature's environmental partners.
    """

    def __init__(self, id, feature_type, atom_numbers, occurrences, envpartners):

        self.id = id
        self.feature_type = feature_type
        self.atom_numbers = atom_numbers
        self.occurrences = occurrences
        self.envpartners = envpartners

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

        return self._data(type="occurrences")

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

        return self._data(type="distances")

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

        superfeature_count = pd.Series({"any": sum(self.occurrences)})
        envpartners_count = pd.Series(
            {envpartner.id: envpartner.count for envpartner in self.envpartners}
        )

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

        return self.count.apply(lambda x: round(x / self.n_frames * 100, 2))

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

        envpartners = pd.DataFrame(
            [getattr(envpartner, type) for envpartner in self.envpartners],
            index=[envpartner.id for envpartner in self.envpartners],
        ).transpose()

        return envpartners
