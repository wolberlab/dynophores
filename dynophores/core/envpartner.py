"""
dynophores.core.envpartner

Handles environmental partner class.
"""


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
    atom_numbers : list of int
        List of atom numbers.
    occurrences : numpy.array
        Occurrences of interaction with environmental partner (0=no, 1=yes) in each frame.
    distances : numpy.array
        Interaction distances in each frame.
    """

    def __init__(
        self, id, residue_name, residue_number, chain, atom_numbers, occurrences, distances
    ):

        if len(occurrences) != len(distances):
            raise ValueError("Occurrences and distances must be of same length.")

        self.id = id
        self.residue_name = residue_name
        self.residue_number = residue_number
        self.chain = chain
        self.atom_numbers = atom_numbers
        self.occurrences = occurrences
        self.distances = distances

    @property
    def n_frames(self):
        """
        Get environmental partner's number of frames.

        Returns
        -------
        int
            Number of frames.
        """

        return len(self.occurrences)

    @property
    def count(self):
        """
        Get number of frames in which the interaction with the environmental partner occurs.

        Returns
        -------
        int
            Count of interaction occurrence.
        """

        return sum(self.occurrences)

    @property
    def frequency(self):
        """
        Get frequency of frames in which the interaction with the environmental partner occurs.

        Returns
        -------
        float
            Frequency of interaction occurrence.
        """

        return round(sum(self.occurrences) / self.n_frames * 100, 2)
