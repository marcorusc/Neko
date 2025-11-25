from __future__ import annotations

from typing import Any, Callable, Iterable, Literal
import pickle
import logging
import pandas as pd
import os
from pypath_common import misc as _common

from ._db.omnipath import omnipath_universe
from ._db import psp as _psp
from ._db import _misc
from ._db import huri as _huri
from ._db import signor as _signor

"""
Access to generic networks from databases, files and standard formats.
"""

_METHODS = {
    'omnipath': omnipath_universe,
}
_REQUIRED_COLS = {
    'source',
    'target',
}

MANDATORY_BOOL_COLS = [
    'is_directed',
    'is_stimulation',
    'is_inhibition',
    'form_complex',
]

MANDATORY_COLUMNS = [
    'source',
    'target',
]
MANDATORY_COLUMNS.extend(MANDATORY_BOOL_COLS)


def network_universe(
        resource: Literal['omnipath'] | pd.DataFrame = 'omnipath',
        **kwargs
    ) -> Universe:
    """
    Generic networks from databases, files and standard formats.

    Args:
        resource:
            Name of the resource or a ready data frame to bypass the built-in
            loading method.
        kwargs:
            Passed to the source specific method. See the specific methods in
            this module for details.

    Note: currently OmniPath PPI is the single available option and serves
    as a placeholder. Later we will dispatch all inputs through this API.
    """

    return Universe(resource, **kwargs)


def omnipath(**kwargs) -> Universe:

    return network_universe('omnipath', **kwargs)


def signor(path: str | None = None, **kwargs) -> Universe:
    # If path is provided and exists, use it. Otherwise, use the new signor() logic to download/process.
    if path and os.path.exists(path):
        return Universe(_signor.signor(path))
    else:
        # Use the new signor() logic: download and process if path is None
        return Universe(_signor.signor())


def phosphosite(
        organism: Literal["human", "mouse", "rat"] = 'human',
        kinase_substrate: str | None = None,
        regulatory_sites: str | None = None,
        **kwargs
    ) -> Universe:

    df = _psp.psp(organism, kinase_substrate, regulatory_sites)

    return Universe(df, name = 'phosphosite')


def huri(dataset: str = 'HI-union') -> Universe:

    df = _huri.huri(dataset)

    return Universe(df, directed = False, name = 'huri')


class Universe:


    def __init__(
            self,
            resources: Literal['omnipath'] | pd.DataFrame = None,
            **param
        ):
        """
        Load and preprocess a generic network from databases or files.
        """

        self._resources = {}
        self._directed = {}
        self._resource = resources
        self.interactions = None
        self.add_resources(resources, **param)
        self.build()


    def add_resources(
            self,
            resources,
            directed: bool = True,
            **param
        ) -> None:
        """
        Add resources to the Universe.

        Parameters:
            resources:
                Can be one of the following:
                - str: Path to a .pickle or .tsv file, or a known resource name
                  like 'omnipath'.
                - pd.DataFrame: A DataFrame containing interaction data.
                - dict: A dictionary mapping resource names to DataFrames.
                - Iterable: A list/tuple of DataFrames.
                - None: No resources to add.
            directed (bool):
                Whether the interactions are directed. Default is True.
            **param:
                Additional parameters:
                - name (str): Name for the resource. Default is '_default'.
                - columns (dict): Column name mapping for renaming columns.

        Returns:
            None
        """
        # Handle None resources gracefully
        if resources is None:
            return

        if isinstance(resources, str):

            if resources.endswith('.pickle'):

                if not os.path.exists(resources):
                    raise FileNotFoundError(
                        f"Pickle file not found: {resources}"
                    )

                with open(resources, 'rb') as fin:

                    resources = pickle.load(fin)

            elif resources.endswith('.tsv'):

                if not os.path.exists(resources):
                    raise FileNotFoundError(f"TSV file not found: {resources}")

                resources = pd.read_csv(resources, sep='\t')

        name = param.get('name', '_default')
        columns = param.get('columns', {})

        if isinstance(resources, str) and resources in _METHODS:

            name = name or resources
            resources = _METHODS[resources](**param)

        if isinstance(resources, pd.DataFrame):

            if resources.empty:
                logging.warning(
                    "Adding empty DataFrame as resource '%s'", name
                )

            self._resources[name] = self._check_columns(
                resources,
                columns,
                directed = directed,
            )
            self._directed[name] = directed

        elif isinstance(resources, dict):

            resources = {
                k: self._check_columns(v, columns)
                for k, v in resources.items()
            }

            if isinstance(directed, bool):

                directed = {k: directed for k in resources.keys()}

            for res_name, resource in resources.items():

                self.add_resources(
                    resource,
                    directed = directed[res_name],
                    name = res_name,
                    **param,
                )

        elif isinstance(resources, Iterable) and not isinstance(resources, str):

            resources_list = list(resources)
            names = [f'{name}_{i}' for i in range(len(resources_list))]

            if isinstance(directed, bool):

                directed = [directed] * len(resources_list)

            for r, d, n in zip(resources_list, directed, names):

                self.add_resources(r, directed = d, name = n, **param)

    def remove_resources(self, name: str) -> bool:
        """
        Remove a resource from the Universe by name.

        Parameters:
            name (str): The name of the resource to remove.

        Returns:
            bool: True if the resource was removed, False if it was not found.
        """
        if name in self._resources:
            del self._resources[name]
            self._directed.pop(name, None)
            # Rebuild the interactions after removing the resource
            self.build()
            return True
        else:
            logging.warning("Resource '%s' not found in Universe", name)
            return False

    def list_resources(self) -> list[str]:
        """
        List all resource names currently loaded in the Universe.

        Returns:
            list[str]: A list of resource names.
        """
        return list(self._resources.keys())

    def get_resource(self, name: str) -> pd.DataFrame | None:
        """
        Get a specific resource DataFrame by name.

        Parameters:
            name (str): The name of the resource.

        Returns:
            pd.DataFrame | None: The resource DataFrame, or None if not found.
        """
        return self._resources.get(name)

    def clear_resources(self) -> None:
        """
        Remove all resources from the Universe.
        """
        self._resources.clear()
        self._directed.clear()
        self.interactions = None

    @staticmethod
    def _check_columns(
            df: pd.DataFrame,
            columns: dict,
            directed: bool = True,
        ) -> pd.DataFrame:

        # If columns is provided, rename the columns of the incoming df
        if columns:
            df = df.rename(columns=columns)

        if 'effect' in df.columns:

            df = _misc.split_effect(df)

        for col in MANDATORY_BOOL_COLS:

            df = _misc.bool_col(df, col)

        # Check if the df contains the required columns
        missing_columns = set(MANDATORY_COLUMNS) - set(df.columns)

        if missing_columns:

            logging.warning("The incoming df is missing some required columns: %s", missing_columns)
            logging.warning("This might lead to issues in running the package.")

        if not directed:

            df = _misc.undirected_to_mutual(df)

        return df


    @staticmethod
    def merge(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
        """
        This function concatenates the provided df with the existing one in the resources object,
        aligning columns and filling in missing data with NaN.

        Parameters:
            df1 (pd.DataFrame): The DataFrame to be added.
            df2 (pd.DataFrame): The DataFrame to be added.

        Raises:
            ValueError: If the 'df' parameter is not a pandas DataFrame.

        Returns:
            None
        """

        # Align columns of both dataframes, filling missing columns with NaN
        all_columns = set(df1.columns).union(set(df2.columns))
        df1 = df1.reindex(columns=all_columns, fill_value=None)
        df2 = df2.reindex(columns=all_columns, fill_value=None)
        df1 = pd.concat([df1, df2])

        return df1.copy()


    def build(self, resources: Iterable[str] | None = None) -> None:

        resources = _common.to_list(resources) or self._resources.keys()

        self.interactions = None

        if not resources:

            return

        for res in resources:
            df = self._resources[res]
            self.interactions = (
                df.copy()
                if self.interactions is None else
                self.merge(self.interactions, df)
            )

        self.interactions.reset_index(drop=True, inplace=True)


    @property
    def resource(self):

        return self._resource if isinstance(self._resource, str) else 'user'


    @property
    def network(self) -> pd.DataFrame:
        """
        The network as it's been read from the original source.
        """

        if not hasattr(self, '_network'):

            self.load()

        return self._network


    @network.setter
    def network(self, value: Any):

       raise AttributeError('The attribute `Universe.network` is read-only.')


    def load(self) -> None:
        """
        Acquire the input data according to parameters.
        """

        self._network = self.method(**self.param)


    @property
    def method(self) -> Callable:
        """
        The method that loads the data.

        ``param`` are to be passed to this method.
        """

        return (
            lambda **kwargs: self._resource
                if isinstance(self._resource, pd.DataFrame) else
            _METHODS[self.resource]
        )


    def __repr__(self) -> str:

        return f'Universe; resources: {", ".join(self._resources.keys())}; size: {len(self)}'


    def __len__(self) -> int:

        return 0 if self.interactions is None else len(self.interactions)


    def check(self) -> bool:
        """
        The network is loaded and contains the mandatory variables.
        """

        return (
            hasattr(self, '_network') and
            not _REQUIRED_COLS - set(self._network.columns)
        )

    @property
    def nodes(self) -> set[str]:

        return (
            set(self.interactions['source']) |
            set(self.interactions['target'])
        )


    def __contains__(self, other) -> bool:

        return other in self.nodes


    def __and__(self, other: set) -> set:

        return self.nodes & other
