# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

from abc import ABC, abstractmethod
from typing import IO, List, Optional, Union

from treezy.tree import Tree


class TreeReader(ABC):
    """Base class for reading trees from a file or stream.

    This class provides a common interface for reading trees in various formats.
    It can read raw strings from the input stream, which can be used for comments
    or other non-tree content.
    """

    def __init__(
        self, path_or_stream: Union[str, IO], taxon_names: Optional[List[str]] = None
    ):
        """Initializes the TreeReader with a file path or an IO stream.

        Parameters
        ----------
        path_or_stream
            a string representing the file path or an IO stream
            (like StringIO or BytesIO) to read trees from.
        taxon_names
            an optional list of taxon names to use for parsing trees.

        Raises
        ------
        TypeError
            if path_or_stream is neither a string nor an IO stream.
        """
        if isinstance(path_or_stream, str):
            self.in_stream = open(path_or_stream, 'r')
            self._path = path_or_stream
        # duck-typing to check if it has a 'readline' method
        elif hasattr(path_or_stream, 'readline'):
            self.in_stream = path_or_stream
            self._path = None
        else:
            raise TypeError("path_or_stream must be a string or an IO stream.")
        self.taxon_names = taxon_names if taxon_names is not None else []
        self._count = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._path and self.in_stream:
            self.in_stream.close()

    def close(self) -> None:
        """Closes the input stream if it was opened by this reader.

        Do not use this method if the input stream was provided externally
        or using a context manager.
        """
        if self.in_stream and not self._path:
            self.in_stream.close()

    @abstractmethod
    def count_trees(self) -> int:
        """Counts the number of trees in the input stream.

        This method reads through the input stream to count the number of trees
        without parsing them.

        Returns
        -------
        int
            the number of trees in the input stream.
        """

    @abstractmethod
    def next(self) -> Optional[Tree]:
        """Returns the next tree from the input stream.

        If there are no more trees, returns None.
        """

    @abstractmethod
    def has_next(self) -> bool:
        """Checks if there is another tree available in the input stream.

        This method reads through the input stream to determine if there is
        another tree available without consuming it.

        Returns
        -------
        bool
            True if there is another tree available, False otherwise.
        """

    @abstractmethod
    def skip_next(self) -> None:
        """Skips the next tree in the input stream.

        This method discards the next tree without parsing it.
        If there are no more trees, it does nothing.
        """

    def parse(self) -> List[Tree]:
        """Parses all trees from the input stream.

        This method reads through the input stream and returns a list of
        Tree objects representing the trees found in the input stream.

        Returns
        -------
        List[Tree]
            a list of Tree objects parsed from the input stream.
        """
        trees = []
        while self.has_next():
            tree = self.next()
            if tree:
                trees.append(tree)
        return trees


class TreeWriter(ABC):
    """Base class for writing trees to a file or stream.

    This class provides a common interface for writing trees in various formats.
    It can write raw strings to the output stream, which can be used for comments
    or other non-tree content.
    """

    def __init__(self, path_or_stream: Union[str, IO], mode: str = 'w'):
        if isinstance(path_or_stream, str):
            self._out_stream = open(path_or_stream, mode)
            self._path = path_or_stream
        # duck-typing to check if it has a 'write' method
        elif hasattr(path_or_stream, 'write'):
            self._out_stream = path_or_stream
            self._path = None
        else:
            raise TypeError("path_or_stream must be a string or an IO stream.")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._path and self._out_stream:
            self._out_stream.close()

    def close(self) -> None:
        """Closes the output stream if it was opened by this writer.

        Do not use this method if the output stream was provided externally
        or using a context manager.
        """
        if self._out_stream and not self._path:
            self._out_stream.close()

    def write_string(self, string: str) -> None:
        """Writes a string to the output stream.

        This method is used to write raw strings directly to the tree file.
        It can be used for comments or other non-tree content.

        Parameters
        ----------
        string
            the string to write to the output stream.
        """
        self._out_stream.write(string + "\n")

    @abstractmethod
    def write(self, trees: Union[Tree, List[Tree]]):
        """Writes one or more trees to the output stream.

        Parameters
        ----------
        trees
            a single Tree object or a list of Tree objects to write.
        """

    @classmethod
    @abstractmethod
    def save(
        cls, path_or_stream: Union[str, IO], trees: Union[Tree, List[Tree]], **options
    ):
        """Saves one or more trees to a file or stream.

        Parameters
        ----------
        path_or_stream
            the file path or IO stream to save the trees to.
        trees
            a single Tree object or a list of Tree objects to save.
        options
            additional options for saving the trees, such as formatting or
                annotations.
        """
