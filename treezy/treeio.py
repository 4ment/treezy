# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

from abc import ABC, abstractmethod
from typing import IO, List, Optional, Union

from treezy.tree import Tree


class TreeReader(ABC):
    def __init__(self, in_stream: IO[str], taxon_names: Optional[List[str]] = None):
        self.in_stream = in_stream
        self.taxon_names = taxon_names if taxon_names is not None else []
        self._count = 0

    @abstractmethod
    def count_trees(self) -> int:
        pass

    @abstractmethod
    def next(self) -> Optional[Tree]:
        pass

    @abstractmethod
    def has_next(self) -> bool:
        pass

    @abstractmethod
    def skip_next(self) -> None:
        pass

    def parse(self) -> List[Tree]:
        trees = []
        while self.has_next():
            tree = self.next()
            if tree:
                trees.append(tree)
        return trees


class TreeWriter:
    def __init__(self, path):
        self._path = path

    def write(self, trees):
        raise NotImplementedError("Subclasses must implement `write`.")

    @classmethod
    def save(cls, path: Union[str, IO], trees: Union[Tree, List[Tree]], **options):
        raise NotImplementedError("Subclasses must implement `save`.")
