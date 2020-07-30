"""
Contains base functions/classes.
"""

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

FEATURE_COLORS = {'HBA': 'firebrick', 'HBD': 'green', 'H': 'gold', 'AR': 'mediumblue', 'PI': 'blue', 'NI': 'red'}
plt.style.use('seaborn')


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
            index=[superfeature.id for superfeature in self.superfeatures]
        ).transpose()

        return occurrence_superfeatures

    @property
    def envpartners_occurrences(self):
        """
        For each superfeature, get its environmental partners' occurrences per environmental partner and frame.

        Returns
        -------
        dict of pandas.DataFrame
            For each superfeature, occurrences (0=no, 1=yes) of an environmental partner (columns) in each frame (row).
        """

        return self._get_envpartners_data(type='occurrences')

    @property
    def envpartners_distances(self):
        """
        For each superfeature, get its environmental partners' distances per environmental partner and frame.

        Returns
        -------
        dict of pandas.DataFrame
            For each superfeature, distances to an environmental partner (columns) in each frame (row).

        """

        return self._get_envpartners_data(type='distances')

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
        Get number of frames in which each dynophore occurs, including the superfeatures and superfeatures'
        environmental partners occurrences.

        Returns
        -------
        pandas.DataFrame
            Dynophore count: The DataFrame shows interaction (yes/no) for superfeatures (rows) to each single
            environmental partner as well as any environmental partner (columns).
        """

        dynophore_count = pd.DataFrame(
            {superfeature.id: superfeature.count for superfeature in self.superfeatures}
        )
        dynophore_count.fillna(0, inplace=True)
        dynophore_count = dynophore_count.astype('int32')

        return dynophore_count

    @property
    def frequency(self):
        """
        Get frequency of frames in which each dynophore occurs, including the superfeatures and superfeatures'
        environmental partners occurrences.

        Returns
        -------
        pandas.DataFrame
            Dynophore frequency: The DataFrame shows interaction (yes/no) for superfeatures (rows) to each single
            environmental partner as well as any environmental partner (columns).
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
            raise KeyError(f'Superfeature name {superfeature_name} is unknown.')

    def from_file(self, dynophore_path):
        """
        Load dynophore data from file to Dynophore instance.

        Parameters
        ----------
        dynophore_path : pathlib.Path
            Path to folder with dynophore data.
        """

        # Get all files and filename components
        dynophore_files = [file for file in dynophore_path.glob('*')]
        dynophore_files_components = [self._get_file_components(i) for i in dynophore_files]

        # Iterate over superfeatures
        superfeatures = []
        for superfeature_file_components in [i for i in dynophore_files_components if i['envpartner_id'] is None]:

            # Iterate over environmental partners
            envpartners = []
            for envpartner_file_components in [
                i for i in dynophore_files_components
                if i['superfeature_id'] == superfeature_file_components['superfeature_id'] and
                i['envpartner_id'] is not None
            ]:
                # Set environmental partner
                envpartner = EnvPartner(
                    envpartner_file_components['envpartner_id'],
                    envpartner_file_components['envpartner_residue_name'],
                    envpartner_file_components['envpartner_residue_number'],
                    envpartner_file_components['envpartner_chain'],
                    envpartner_file_components['envpartner_atom_numbers'],
                    np.loadtxt(fname=envpartner_file_components['filepath'], dtype=int, delimiter=',', usecols=1),
                    np.loadtxt(fname=envpartner_file_components['filepath'], dtype=float, delimiter=',', usecols=0)
                )
                envpartners.append(envpartner)

            # Set superfeature
            superfeature = Superfeature(
                superfeature_file_components['superfeature_id'],
                superfeature_file_components['superfeature_feature_type'],
                superfeature_file_components['superfeature_atom_numbers'],
                np.loadtxt(fname=superfeature_file_components['filepath'], dtype=int),
                envpartners
            )
            superfeatures.append(superfeature)

        # Add superfeature to dynophore
        self.id = dynophore_files_components[0]['dynophore_id']
        self.superfeatures = superfeatures

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

    def plot_superfeatures_occurrences(self, superfeature_names=None, color_by_feature_type=True, max_frames=1000):
        """
        Plot the superfeatures' occurrences as barcode.

        Parameters
        ----------
        superfeature_names : None or str or list of str
            Select superfeatures to plot or select all (default).
        color_by_feature_type : bool
            Color barcode by feature type (default) or color all in black.
        max_frames : int
            Number of frames to display in barcode plot. If input data contains more than `max_frames`,
            `max_frames` equidistant frames will be selected.
        """

        occurrences = self.superfeatures_occurrences

        if superfeature_names is not None:
            if isinstance(superfeature_names, str):
                superfeature_names = [superfeature_names]

            # Get all superfeature names that are in data
            superfeature_names_curated = [i for i in superfeature_names if i in occurrences.columns]
            superfeature_names_omitted = list(set(superfeature_names)-set(superfeature_names_curated))
            if len(superfeature_names_omitted) > 0:
                print(f'Superfeature names {superfeature_names_omitted} omitted because unknown.')

            # Select subset
            occurrences = occurrences[superfeature_names_curated]

        # Prepare data
        occurrences = self._prepare_plot_occurrences(occurrences, max_frames)

        # Feature type colors?
        if color_by_feature_type:
            feature_types = [i.split('[')[0] for i in occurrences.columns]
            colors = [FEATURE_COLORS[i] if i in FEATURE_COLORS.keys() else 'black' for i in feature_types]
        else:
            colors = None

        # Plot (plot size depending on number barcodes)
        fig, ax = plt.subplots(figsize=(10, occurrences.shape[1] / 2))
        occurrences.plot(marker='.', markersize=5, linestyle='', legend=None, ax=ax, color=colors)
        # Set y tick labels
        ax.set_yticks(range(0, occurrences.shape[1]+2))
        ax.set_yticklabels([''] + occurrences.columns.to_list() + [''])
        ax.invert_yaxis()
        # Set x axis limits and label
        ax.set_xlabel('frame index')
        ax.set_xlim((occurrences.index[0], occurrences.index[-1]))

    def plot_envpartners_occurrences(self, superfeature_name, max_frames=1000):
        """
        Plot a superfeature's interaction ocurrences with its interaction partners.

        Parameters
        ----------
        superfeature_name : str
            Superfeature name
        max_frames : int
            Number of frames to display in barcode plot. If input data contains more than `max_frames`,
            `max_frames` equidistant frames will be selected.
        """

        self.is_superfeature(superfeature_name)

        # Prepare data
        occurrences = self._prepare_plot_occurrences(self.envpartners_occurrences[superfeature_name], max_frames)

        # Plot (plot size depending on number barcodes)
        fig, ax = plt.subplots(figsize=(10, occurrences.shape[1] / 2))
        occurrences.plot(marker='.', markersize=5, linestyle='', legend=None, ax=ax)
        # Set y tick labels
        ax.set_yticks(range(0, occurrences.shape[1]+2))
        ax.set_yticklabels([''] + occurrences.columns.to_list() + [''])
        ax.invert_yaxis()
        # Set x axis limits and label
        ax.set_xlabel('Frame')
        ax.set_xlim((occurrences.index[0], occurrences.index[-1]))

    def plot_envpartners(self, superfeature_name, max_frames=1000):
        """
        Plot interaction data for a superfeature, i.e. occurrences (frame series) and distances (frame series and
        histogram).

        Parameters
        ----------
        superfeature_name : str
            Superfeature name
        max_frames : int
            Number of frames to display in barcode plot. If input data contains more than `max_frames`,
            `max_frames` equidistant frames will be selected.
        """

        self.is_superfeature(superfeature_name)
        occurrences = self._prepare_plot_envparters_occurrences(superfeature_name, max_frames)
        distances = self._prepare_plot_envpartners_distances(superfeature_name, max_frames)

        # Set up plot
        fig, axes = plt.subplots(
            nrows=2,
            ncols=2,
            figsize=(15, 7),
            sharey='row',
            sharex='col',
            gridspec_kw={
                'width_ratios': [3, 1],
                'wspace': 0.05,
                'hspace': 0.05
            }
        )

        # Subplot (0, 0): Interaction occurrences (barplot)
        occurrences.plot(ax=axes[0][0], kind='line', legend=None, marker='.', linestyle='')
        # Set y tick labels (tick per envpartner but do not show label)
        axes[0][0].set_yticks(range(0, occurrences.shape[1]+2))
        axes[0][0].set_yticklabels('')
        # Set x axis limits and label
        axes[0][0].set_xlabel('frame index')
        axes[0][0].set_xlim((occurrences.index[0], occurrences.index[-1]))
        # Show interactions from top to bottom from most common to rarest interactions
        axes[0][0].invert_yaxis()

        # Subplot (0, 1): Empty (will hold legend from subplot (1, 1))
        axes[0][1].axis('off')

        # Subplot (1, 0): Distance time series
        distances.plot(ax=axes[1][0], kind='line', legend=None, linewidth=1)
        axes[1][0].set_xlabel('Frame index', fontsize=16)
        axes[1][0].set_ylabel(r'Distance [$\AA$]', fontsize=16)

        # Subplot (1, 1): Distance histogram
        bins = range(0, math.ceil(distances.max().max())+1)
        distances.plot(
            ax=axes[1][1],
            kind='hist',
            bins=bins,
            orientation="horizontal",
            alpha=0.7,
            density=True,
            xlim=(0, 1)
        )
        axes[1][1].set_xlabel('Frequency', fontsize=16)
        axes[1][1].legend(loc=6, bbox_to_anchor=(0, 1.5), fontsize=12)

    def plot_superfeatures_vs_envpartners(self):
        """
        Plot heatmap of interactions between superfeatures and interaction partners.
        """

        # Sort superfeatures by overall frequency
        data = self.frequency[self.frequency.loc['any', :].sort_values(ascending=False).index]
        data.fillna(0, inplace=True)
        sns.heatmap(data, cmap='Blues')

    def plot_envpartner_distances(self, superfeature_name, kind):
        """
        Plot interaction distances for a superfeatures as frame series or histogram.

        Parameters
        ----------
        superfeature_name : str
            Superfeature name.
        kind : str
            Plot kind, 'line' (distance vs. frames) or 'hist' (distance histogram)
        """

        data = self.envpartners_distances[superfeature_name]

        if kind == 'line':
            fig, ax = plt.subplots(figsize=(10, 5))
            data.plot(kind='line', ax=ax)
            ax.set_xlim((0, data.shape[0]))
            ax.set_xlabel('Frame index')
            ax.set_ylabel(r'Distance [$\AA$]')
        elif kind == 'hist':
            ax = data.plot(kind='hist')
            ax.set_xlabel(r'Distance [$\AA$]')
        else:
            raise ValueError('Plotting kind is unknown. Choose from "line" and "hist".')

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

    def _get_envpartners_data(self, type='occurrences'):
        """
        Get occurrences or distances of all superfeatures' environmental partners.

        Parameters
        ----------
        type : str
            Data type: occurrences (default) or distances.

        Returns
        -------
        dict of DataFrame
            Occurrences (default) or distances for a superfeature's (dict key) environmental partners (columns) for all
            frames (rows).
        """

        types = ['occurrences', 'distances']
        if type in types:
            pass
        else:
            raise ValueError(f'Wrong type. Select from: {", ".join(types)}')

        envpartners = {}

        for superfeature in self.superfeatures:
            superfeature_envpartners = pd.DataFrame(
                [getattr(envpartner, type) for envpartner in superfeature.envpartners],
                index=[envpartner.id for envpartner in superfeature.envpartners]
            ).transpose()

            envpartners[superfeature.id] = superfeature_envpartners

        return envpartners

    @staticmethod
    def _prepare_plot_occurrences(occurrences, max_frames=1000):
        """
        Prepare data for plotting occurrences (superfeatures or superfeature interactions):
        - Sort by interaction frequency
        - Get subset of data if more than 1000 frames
        - Transform 1 in binary values to rank in plot

        Parameters
        ----------
        occurrences : pandas.DataFrame
            Occurrences (0 or 1) per frame (rows) for one or more superfeatures or superfeature interactions (columns).
        max_frames : int
            Number of frames to display in barcode plot. If input data contains more than `max_frames`,
            `max_frames` equidistant frames will be selected.

        Returns
        -------
        pandas.DataFrame
            Occurrences, ready for plotting.
        """

        # Sort data by ratio
        ratio = round(occurrences.apply(sum) / occurrences.shape[0] * 100, 2).sort_values(ascending=False)
        occurrences = occurrences[ratio.index]

        # Get subset of data if more than 1000 frames
        if occurrences.shape[0] > 1000:
            selected_indices = [i for i in range(0, 1000, math.floor(1002/max_frames))]
            occurrences = occurrences.iloc[selected_indices, :]

        # Transform 1 in binary values to rank in plot
        occurrences_plot = {}
        for i, (name, data) in enumerate(occurrences.iteritems()):
            data = data.replace([0, 1], [None, i+1])
            occurrences_plot[name] = data
        occurrences_plot = pd.DataFrame(occurrences_plot)

        return occurrences_plot

    def _prepare_plot_envparters_occurrences(self, superfeature_name, max_frames=1000):
        """
        Prepare data for plotting superfeature interaction occurrences.

        Parameters
        ----------
        superfeature_name : str
            Superfeature name
        max_frames : int
            Number of frames to display in barcode plot. If input data contains more than `max_frames`,
            `max_frames` equidistant frames will be selected.

        Returns
        -------
        pandas.DataFrame
            Occurrences, ready for plotting.
        """

        occurrences = self.envpartners_occurrences[superfeature_name]
        occurrences = self._prepare_plot_occurrences(occurrences, max_frames)

        return occurrences

    def _prepare_plot_envpartners_distances(self, superfeature_name, max_frames=1000):
        """
        Prepare data for plotting interaction distances for a superfeature:
        - Sort by interaction frequency

        Parameters
        ----------
        superfeature_name : str
            Superfeature name
        max_frames : int
            Number of frames to display in barcode plot. If input data contains more than `max_frames`,
            `max_frames` equidistant frames will be selected.

        Returns
        -------
        pandas.DataFrame
            Interaction distances for a superfeature, ready for plotting.
        """

        # Get data
        distances = self.envpartners_distances[superfeature_name]

        # Sort data by ratio
        ratio = self.frequency.loc[:, superfeature_name].dropna().drop('any').sort_values(ascending=False)
        distances = distances[ratio.index]

        return distances


class Superfeature:
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
    occurrences : pandas.DataFrame
        Per superfeature (columns), occurrences (0=no, 1=yes) in each frame (row).
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
        Get the superfeature's environmental partners' occurrences per environmental partner and frame.

        Returns
        -------
        pandas.DataFrame
            Occurrences (0=no, 1=yes) of an environmental partner (columns) in each frame (row).

        """

        return self._get_envpartners_data(type='occurrences')

    @property
    def envpartners_distances(self):
        """
        Get the superfeature's environmental partners' distances per environmental partner and frame.

        Returns
        -------
        pandas.DataFrame
            Distances to an environmental partner (columns) in each frame (row)
        """

        return self._get_envpartners_data(type='distances')

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
        Get number of frames in which the superfeature occurs, including the superfeature's environmental partners
        occurrences.

        Returns
        -------
        pandas.Series
            Superfeature count: The Series shows interactions (yes/no) to each single environmental partner as well as
            any environmental partner.
        """

        superfeature_count = pd.Series({'any': sum(self.occurrences)})
        envpartners_count = pd.Series(
            {envpartner.id: envpartner.count for envpartner in self.envpartners}
        )

        return superfeature_count.append(envpartners_count)

    @property
    def frequency(self):
        """
        Get frequency of frames in which the superfeature occurs, including the superfeature's environmental partners
        occurrences.

        Returns
        -------
        pandas.Series
            Superfeature frequency: The Series shows interactions (yes/no) to each single environmental partner as well
            as any environmental partner.
        """

        return self.count.apply(lambda x: round(x / self.n_frames * 100, 2))

    def _get_envpartners_data(self, type='occurrences'):
        """
        Get occurrences or distances of a superfeature's environmental partners.

        Parameters
        ----------
        type : str
            Data type: occurrences (default) or distances.

        Returns
        -------
        DataFrame
            Occurrences (default) or distances of a superfeature's environmental partners (columns) for all
            frames (rows).
        """

        types = ['occurrences', 'distances']
        if type in types:
            pass
        else:
            raise ValueError(f'Wrong type. Select from: {", ".join(types)}')

        envpartners = pd.DataFrame(
            [getattr(envpartner, type) for envpartner in self.envpartners],
            index=[envpartner.id for envpartner in self.envpartners]
        ).transpose()

        return envpartners


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

    def __init__(self, id, residue_name, residue_number, chain, atom_numbers, occurrences, distances):

        if len(occurrences) != len(distances):
            raise ValueError('Occurrences and distances must be of same length.')

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


class Plotting:
    """
    Class to store plotting functions
    """

    def __init__(self):
        pass

    @staticmethod
    def plot_occurrences(occurrences, color_by_feature_type=True, max_frames=1000):
        """
        Create a barcode plot for superfeature or interaction occurrences.

        Parameters
        ----------
        occurrences : pandas.DataFrame
            Occurrences (0 or 1) per frame (rows) for one or more superfeatures or interactions (columns).
        color_by_feature_type : bool
            Color barcode by feature type (default) or color all in black.
        max_frames : int
            Number of frames to display in barcode plot. If input data contains more than `max_frames`,
            `max_frames` equidistant frames will be selected.
        """
        # Set plotting style
        plt.style.use('seaborn-dark')

        # Sort data by ratio
        ratio = round(occurrences.apply(sum) / occurrences.shape[0] * 100, 2)
        ratio.sort_values(inplace=True)
        occurrences = occurrences[ratio.index]

        # Get subset of data if more than 1000 frames
        if occurrences.shape[0] > 1000:
            selected_indices = [i for i in range(0, 1000, math.floor(1002/max_frames))]
            occurrences = occurrences.iloc[selected_indices, :]

        # Transform 1 in binary values to rank in plot
        occurrences_plot = {}
        for i, (name, data) in enumerate(occurrences.iteritems()):
            data = data.replace([0, 1], [None, i+1])
            occurrences_plot[name] = data
        occurrences_plot = pd.DataFrame(occurrences_plot)

        # Feature type colors?
        if color_by_feature_type:
            feature_types = [i.split('[')[0] for i in occurrences_plot.columns]
            colors = [FEATURE_COLORS[i] if i in FEATURE_COLORS.keys() else 'black' for i in feature_types]
        else:
            colors = None

        # Plot (plot size depending on number barcodes)
        fig, ax = plt.subplots(figsize=(10, occurrences_plot.shape[1] / 2))
        occurrences_plot.plot(marker='|', markersize=5, linestyle='', legend=None, ax=ax, color=colors)
        # Set y tick labels
        ax.set_yticks(range(0, occurrences_plot.shape[1]+2))
        ax.set_yticklabels([''] + occurrences_plot.columns.to_list() + [''])
        # Set x axis limits and label
        ax.set_xlabel('frame index')
        ax.set_xlim((occurrences_plot.index[0], occurrences_plot.index[-1]))
