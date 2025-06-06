#!/usr/bin/env python

#
# This file is part of the `neko` Python module
#
# Copyright 2023
# Insitut Curie
#
# File author(s): Marco Ruscone (marco.ruscone@curie.fr)
#
# Distributed under the GPLv3 license
# See the file `LICENSE` or read a copy at
# https://www.gnu.org/licenses/gpl-3.0.txt
#

"""
Package metadata (version, authors, etc).
"""

__all__ = ['get_metadata']

import os
import pathlib
import importlib.metadata

import toml

_VERSION = '0.0.1'


def get_metadata():
    """
    Basic package metadata.

    Retrieves package metadata from the current project directory or from
    the installed package.
    """

    here = pathlib.Path(__file__).parent
    pyproj_toml = 'pyproject.toml'
    meta = {}

    for project_dir in (here, here.parent):

        toml_path = str(project_dir.joinpath(pyproj_toml).absolute())

        if os.path.exists(toml_path):

            pyproject = toml.load(toml_path)

            meta = {
                'name': pyproject['tool']['poetry']['name'],
                'version': pyproject['tool']['poetry']['version'],
                'author': pyproject['tool']['poetry']['authors'],
                'license': pyproject['tool']['poetry']['license'],
                'full_metadata': pyproject,
            }

            break

    if not meta:

        try:

            meta = {
                k.lower(): v for k, v in
                importlib.metadata.metadata(here.name).items()
            }

        except importlib.metadata.PackageNotFoundError:

            pass

    meta['version'] = meta.get('version', None) or _VERSION

    return meta


metadata = get_metadata()
__version__ = metadata.get('version', None)
__author__ = ["Marco Ruscone", "Eirini Tsirvouli", "Andrea Checcoli", "Dénes Turei", "Aasmund Flobak", "Emmanuel Barillot", "Loredana Martignetti", "Julio Saez-Rodriguez", "Laurence Calzone"]
__license__ = metadata.get('license', None)
