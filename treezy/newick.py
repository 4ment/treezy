# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

from typing import IO, List, Optional

from treezy.tree import Tree
from treezy.treeio import TreeReader


class NewickReader(TreeReader):
    def __init__(self, in_stream: IO[str], taxon_names: Optional[List[str]] = None):
        super().__init__(in_stream, taxon_names)
        self._buffer = None

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
            tree = Tree.from_newick(self._buffer, self.taxon_names)
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
                    else:
                        self._buffer = buffer.strip()
                elif buffer == '':
                    self._buffer = None
                    break
