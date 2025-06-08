# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

from typing import IO, List, Optional, Union

from treezy.tree import Tree
from treezy.treeio import TreeReader


class NewickReader(TreeReader):
    """A class for reading trees in Newick format from a file or stream.

    This class extends `TreeReader` to handle Newick-formatted trees.
    Each tree is expected to be represented in standard Newick notation.

    A typical file might contain::

    ((A,B),C);
    ((A,C),B);

    Example
    -------
    >>> import io
    >>> f = io.StringIO("((A,B),C);\\n((A,C),B);")
    >>> with NewickReader(f) as reader:
    ...     for tree in reader.parse():
    ...         print(tree.newick())
    ((A,B),C);
    ((A,C),B);
    """

    def __init__(
        self,
        path_or_stream: Union[str, IO],
        taxon_names: Optional[List[str]] = None,
        **options,
    ):
        """Initializes the NewickReader with a file path or an IO stream.

        Parameters
        ----------
        path_or_stream
            a string representing the file path or an IO stream
            (like StringIO or BytesIO) to read trees from.
        taxon_names
            an optional list of taxon names to use for parsing trees.
        options
            Additional keyword arguments. Currently supports:
                - `strip_quotes`: bool, whether to strip quotes from taxon names
                  (default: False)

        Raises
        ------
        TypeError
            if path_or_stream is neither a string nor an IO stream.
        """
        super().__init__(path_or_stream, taxon_names)
        self._buffer = None
        self._options = options

    def count_trees(self) -> int:
        if self._count > 0:
            return self._count

        self._count = 0
        while True:
            line = self.in_stream.readline()
            if line == '':
                break
            if line.startswith("("):
                self._count += 1

        self.in_stream.seek(0)
        return self._count

    def next(self) -> Optional[Tree]:
        if self.has_next():
            tree = Tree.from_newick(self._buffer, self.taxon_names, **self._options)
            assert isinstance(
                tree, Tree
            ), f"Parsed object is not a Tree instance {self._buffer}"
            self._buffer = None
            return tree
        else:
            return None

    def has_next(self) -> bool:
        if self._buffer is not None:
            return True

        while True:
            self._buffer = self.in_stream.readline()
            if self._buffer == '':
                self._buffer = None
                break
            if self._buffer.startswith("("):
                self._buffer = self._buffer.strip()
                break

        return self._buffer is not None

    def skip_next(self) -> None:
        # discard the tree in the buffer that has been used
        if self._buffer is not None:
            self._buffer = None
        else:
            # discard the first tree
            while True:
                buffer = self.in_stream.readline()
                if buffer.startswith("("):
                    if self._buffer is not None:
                        break
                    self._buffer = buffer.strip()
                elif buffer == '':
                    self._buffer = None
                    break
