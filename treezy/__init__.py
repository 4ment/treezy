# Copyright 2026 Mathieu Fourment
# SPDX-License-Identifier: MIT

"""Simple library for manipulating phylogenetic trees.

This package provides tools for working with phylogenetic trees, including
parsing, manipulation, and comparison operations.
"""

__version__ = "0.0.1"

from treezy.newick import NewickReader
from treezy.nexus import NexusReader, NexusWriter
from treezy.node import Node
from treezy.tree import Tree

__all__ = [
    "NewickReader",
    "NexusReader",
    "NexusWriter",
    "Node",
    "Tree",
]
