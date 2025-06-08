# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

import re
from typing import IO, List, Optional, Union

from treezy.tree import Tree
from treezy.treeio import TreeReader, TreeWriter


class NexusReader(TreeReader):
    """A reader for Nexus format files, which can read trees and their associated
    metadata. This class extends TreeReader to handle Nexus-specific parsing, including
    translate blocks and tree comments. It can read multiple trees from a Nexus file
    and provides methods to count trees, parse them, and access their properties.
    It supports reading translate blocks to map shorthand taxon names to full names,
    and can handle comments within tree definitions.
    It also provides methods to parse or skip trees one by one, which helps avoid memory
    overload when working with large datasets or when subsampling trees.

    Example
    -------
    A typical file might contain::

        #NEXUS
        Begin trees;
        Translate
            1 A,
            2 B,
            3 C;
        tree tree1 = [&R] ((1,2),3);
        tree tree2 = [&R] ((1,3),2);
        End;

    Example
    -------
    >>> import io
    >>> nexus_data = '''#NEXUS
    ... Begin trees;
    ... Translate
    ...     1 A,
    ...     2 B,
    ...     3 C;
    ... tree tree1 = [&R] ((1,2),3);
    ... tree tree2 = [&R] ((1,3),2);
    ... End;
    ... '''
    >>> f = io.StringIO(nexus_data)
    >>> with NexusReader(f) as reader:
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
        """Initializes a NexusReader instance.

        Parameters
        ----------
        in_stream
            the file path or IO stream to read the Nexus file.
        taxon_names
            list of taxon names to use for the trees.
        options
            Additional keyword arguments. Currently supports:
                - `strip_quotes`: bool, whether to strip quotes from taxon names in the
                translate block (default: False)
        """
        super().__init__(path_or_stream, taxon_names)
        self._taxon_map = {}  # maps taxon names to their indices in taxon_names
        self._translate_map = None  # maps shorthand to full taxon names
        self._current_tree_string = None
        self._translate_parsed = False
        self._strip_quotes = options.get("strip_quotes", False)
        self._options = options

    def _find_block(self, block_name: str) -> bool:
        """Searches for a specific block in the Nexus file.

        This method scans through the input stream to find a block that starts with
        "Begin <block_name>". It returns True if the block is found, otherwise False.
        Mostly used to find "Begin trees;".

        Parameters
        ----------
        block_name
            the name of the block to search for (e.g., "trees").

        Returns
        -------
        bool
            True if the block is found, False otherwise.
        """
        for line in self.in_stream:
            if line.lstrip().lower().startswith(f"begin {block_name}"):
                return True
        return False

    def count_trees(self) -> int:
        """Counts the number of trees in the Nexus file.

        This method scans through the input stream to count lines that start with
        ``tree`` after the ``Begin trees;`` block. It returns the total count of
        such lines.
        The stream is reset to the beginning after counting.

        Returns
        -------
        int
            the number of trees found in the Nexus file.
        """
        if self._count > 0:
            return self._count

        found = self._find_block("trees")
        if not found:
            return 0

        self._count = 0
        for line in self.in_stream:
            start_line = line.lstrip()[:4].lower()
            if start_line == "end;":
                break
            if start_line == "tree":
                self._count += 1
        self.in_stream.seek(0)
        return self._count

    def _parse_translate(self):
        """Parses the 'Translate' block in the Nexus file."""
        end = False
        self._translate_map = {}

        while not end:
            line = self._next_line_uncommented().strip()
            if line == ";":
                break
            if ";" in line:
                line, _ = line.split(";", 1)
                end = True

            line = line.rstrip(",")
            for part in line.split(","):
                part = part.strip()
                shorthand, name = re.split(r'[ \t]+', part, maxsplit=1)
                if self._strip_quotes:
                    name = name.strip('"\'')
                self._translate_map[shorthand] = name

    def _next_line_uncommented(self) -> str:
        line = self.in_stream.readline()
        line = line.strip()

        if "[" not in line:
            return line

        # Remove comments
        while "[" in line:
            start = line.find("[")
            end = line.find("]", start)
            if end != -1:
                line = line[:start] + line[end + 1 :]
            else:
                # Multiline comment
                rest = self.in_stream.readline()
                rest = rest.strip()
                while "]" not in rest and rest:
                    rest = self.in_stream.readline()
                    rest = rest.strip()
                if "]" in rest:
                    end = rest.find("]")
                    line = line[:start] + rest[end + 1 :]
                else:
                    line = line[:start]
        return line

    def _parse_tree_line(self, line: str) -> Tree:
        """Parses a single tree line from the Nexus file.

        This method extracts the Newick string from a line that starts with `tree`,
        handles comments, and translates taxon names if a "Translate" block was present.

        Parameters
        ----------
        line
            the line containing the tree definition.

        Returns
        -------
        Tree
            a Tree object constructed from the Newick string in the line.

        Raises
        ------
        ValueError
            if the line does not contain a valid tree definition or if
            the Newick string is malformed.
        """
        start = line.lower().find("tree")
        if start == -1:
            raise ValueError("Not a tree line")

        n = len(line)
        idx = 0
        bracket_count = 0
        comments = []
        while idx < n:
            if line[idx] == '[':
                if bracket_count == 0:
                    start = idx
                bracket_count += 1
            elif line[idx] == ']':
                bracket_count -= 1
                if bracket_count == 0:
                    comments.append(line[start : idx + 1])
            elif bracket_count == 0 and line[idx] == '(':
                break
            idx += 1

        if '[&U]' in comments:
            comments.remove('[&U]')
        elif '[&R]' in comments:
            comments.remove('[&R]')

        if idx == n:
            raise ValueError("No '(' found in tree line")
        # Find last ';'
        end = line.rfind(";")
        if end == -1:
            raise ValueError("No ';' found in tree line")
        newick = line[idx : end + 1]

        if self._translate_map is not None:
            tree = Tree.from_newick(newick, [])
            temporary_taxon_names = []
            for node in tree.root.postorder():
                if node.is_leaf:
                    # if node.name in self._translate_map:
                    node.name = self._translate_map[node.name]
                    temporary_taxon_names.append(node.name)

                    # taxon_name not empty: either provided or not the first tree
                    if node.name in self._taxon_map:
                        node.id = self._taxon_map[node.name]
                    # first tree and taxon_name is empty
                    else:
                        self._taxon_map[node.name] = node.id
            if len(self._taxon_map) != 0 and len(self.taxon_names) == 0:
                self.taxon_names = [None] * len(self._taxon_map)
                for name, idx in self._taxon_map.items():
                    self.taxon_names[idx] = name

            if set(self.taxon_names) != set(temporary_taxon_names):
                missing = set(temporary_taxon_names) - set(self.taxon_names)
                extra = set(self.taxon_names) - set(temporary_taxon_names)
                if missing:
                    print(f"Missing taxon names: {missing}")
                if extra:
                    print(f"Extra taxon names: {extra}")
                raise ValueError(
                    "NexusFile: Taxon names in tree do not match provided taxon names"
                )
            tree.taxon_names = self.taxon_names
        else:
            tree = Tree.from_newick(newick, self.taxon_names, **self._options)

        if len(comments) > 0:
            tree.comment = comments[0]
        return tree

    def _point_to_first_tree(self):
        """Points the reader to the first tree in the Nexus file."""
        if not self._translate_parsed:
            if len(self.taxon_names) != 0:
                for i, name in enumerate(self.taxon_names):
                    self._taxon_map[name] = i
            found = self._find_block("trees")
            if not found:
                print("No trees block found in Nexus file")
                return

            for line in self.in_stream:
                if line == '':
                    self._current_tree_string = None
                    break
                line = line.lstrip()
                if line[:4].lower() == "end;":
                    break
                if line[:9].lower() == "translate":
                    self._parse_translate()
                elif line[:4].lower() == "tree":
                    self._current_tree_string = line.rstrip()
                    break
        self._translate_parsed = True

    def has_next(self) -> bool:
        if self._current_tree_string is not None:
            return True

        if not self._translate_parsed:
            self._point_to_first_tree()
            return self._current_tree_string is not None

        while True:
            line = self.in_stream.readline()
            if line == '':
                self._current_tree_string = None
                break
            line = line.lstrip()
            if line[:4].lower() == "end;":
                self._current_tree_string = None
                break
            if line[:4].lower() == "tree":
                self._current_tree_string = line.rstrip()
                break
        return self._current_tree_string is not None

    def next(self) -> Optional[Tree]:
        """Returns the next tree from the Nexus file.

        This method reads the next tree line from the Nexus file, parses it into a
        :class:`~treezy.tree.Tree` object, and returns it. If there are no more trees,
        it returns ``None``.

        Returns
        -------
        Optional[Tree]:
            the next tree object, or `None` if no more trees are available.
        """
        if self.has_next():
            tree = self._parse_tree_line(self._current_tree_string)
            self._current_tree_string = None
            return tree
        return None

    def skip_next(self) -> None:
        """Skips the next tree in the Nexus file."""
        if self.has_next():
            while True:
                line = self.in_stream.readline()
                if line == '':
                    self._current_tree_string = None
                    break
                line = line.lstrip()
                if line[:4].lower() == "tree":
                    self._current_tree_string = line.rstrip()
                    break


class NexusWriter(TreeWriter):
    """A writer for Nexus format files, which can write trees and their associated
    metadata. This class extends TreeWriter to handle Nexus-specific formatting,
    including translate blocks and tree comments. It can write multiple trees to a Nexus
    file and provides methods to start and end blocks, write translate blocks, and
    write trees in Newick format. It supports options for formatting, such as including
    a translation block, adding comments for each tree, and specifying a prefix for tree
    names. The writer can be used as a context manager to ensure proper file handling.
    """

    def __init__(self, path_or_stream: Union[str, IO], **options):
        """Initializes a NexusWriter instance.

        Parameters
        ----------
        path
            the file path to save the Nexus file.
        options
            formatting options for the Nexus file.
                - `mode`: str, file mode for writing (default: 'w')
                - `include_translate`: bool, whether to include a translation block
                  (default: True)
                - `include_tree_comment`: bool, whether to include comments for each
                  tree (default: False)
                - `tree_prefix`: str, prefix for tree names (default: `STATE_`)
                - Other Newick options as needed.
        """
        super().__init__(path_or_stream)
        self.options = options
        self._tree_counter = 0
        self._mode = options.get('mode', 'w')
        self._translate_written = False
        nexus_keys = {
            'include_translate',
            'include_tree_comment',
            'mode',
            'tree_prefix',
        }
        self._newick_options = {k: v for k, v in options.items() if k not in nexus_keys}

    def __enter__(self):
        self._out_stream.write("#NEXUS\n")
        return self

    def begin_block(self, block_name: str) -> None:
        """Begins a new block in the Nexus file.

        This method writes the opening line for a block, such as ``Begin trees;`` or
        ``Begin taxa;``.

        Parameters
        ----------
        block_name
            the name of the block to begin (e.g., "trees", "taxa").
        """
        self._out_stream.write(f"Begin {block_name};\n")

    def end_block(self) -> None:
        """Ends the current block in the Nexus file."""
        self._out_stream.write("End;\n")

    @staticmethod
    def write_translate_and_map(writer, taxon_names: List[str]) -> dict:
        """Writes the ``Translate`` block to the Nexus file.

        The translate block maps taxon names to their indices (starting from 1)
        in the taxon_names list.

        Example
        -------
        A taxon_names list of `['B', 'A', 'C']` will result in::

            Translate
                1 B,
                2 A,
                3 C;

        Parameters
        ----------
        taxon_names:
            List of taxon names to include in the translate block.
        """
        writer.write("Translate\n")
        translate_list = []
        translate_taxa = {}
        for idx, taxon in enumerate(taxon_names):
            translate_list.append(f"  {idx+1} {taxon}")
            translate_taxa[taxon] = str(idx + 1)
        writer.write(",\n".join(translate_list) + ";\n")
        return translate_taxa

    def write_taxa_block(self, taxon_names: List[str]) -> None:
        """Writes the ``Taxa`` block to the Nexus file.

        Parameters
        ----------
        taxon_names
            List of taxon names to include in the taxa block.
        """
        self._out_stream.write("Begin taxa;")
        self._out_stream.write(f"  Dimensions ntax={len(taxon_names)};")
        self._out_stream.write("  Taxlabels " + " ".join(taxon_names) + ";")
        self._out_stream.write("End;")

    def write_translate(self, taxon_names: List[str]) -> None:
        """Writes the ``Translate`` block to the Nexus file.

        The translate block maps taxon names to their indices (starting from 1)
        in the taxon_names list.

        Example
        -------
        A taxnon_names list of `['A', 'B']` will result in:

            Translate
                1 A,
                2 B;

        Parameters
        ----------
        taxon_names
            List of taxon names to include in the translate block.
        """
        self._newick_options['translator'] = NexusWriter.write_translate_and_map(
            self._out_stream, taxon_names
        )
        self._translate_written = True

    def write(self, trees: Union[Tree, List[Tree]]) -> None:
        """Writes trees in Nexus format to the output stream.

        This method writes a single tree or a list of trees in Nexus format,
        including a translation block if specified in the options. It handles
        formatting options such as including comments for each tree and
        specifying a prefix for tree names.

        Parameters
        ----------
        trees
            a single Tree object or a list of Tree objects to write.
        """
        if isinstance(trees, str):
            self._out_stream.write(trees + "\n")
        elif not isinstance(trees, list):
            trees = [trees]

        include_translate = self.options.get('include_translate', True)

        if include_translate and not self._translate_written:
            self.write_translate(trees[0].taxon_names)

        include_tree_comment = self.options.get('include_tree_comment', False)
        tree_prefix = self.options.get('tree_prefix', 'STATE_')

        for tree in trees:
            self._tree_counter += 1
            name = f"{tree_prefix}{self._tree_counter}"
            newick = tree.newick(**self._newick_options)
            rooting = 'R' if tree.is_rooted else 'U'
            self._out_stream.write(f"tree {name} ")
            if include_tree_comment and tree.comment:
                self._out_stream.write(f"{tree.comment} ")
            self._out_stream.write(f"= [&{rooting}] {newick}\n")

    @classmethod
    def save(
        cls, path_or_stream: Union[str, IO], trees: Union[Tree, List[Tree]], **options
    ):
        """Saves trees in Nexus format to the specified path.

        The `options` dictionary can include:
            - `mode`: str, file mode for writing (default: 'w')
            - `include_translate`: bool, whether to include a translation block
              (default: ``True``)
            - `include_tree_comment`: bool, whether to include comments for each tree
              (default: ``False``)
            - `tree_prefix`: str, prefix for tree names (default: ``STATE_``)
            - Other Newick options as needed.

        Parameters
        ----------
        path_or_stream
            the file path or IO stream to save the Nexus file.
        trees
            a single Tree object or a list of tree objects to save.
        options
            formatting options for the Nexus file.
        """
        mode = options.get('mode', 'w')
        include_translate = options.get('include_translate', True)
        include_tree_comment = options.get('include_tree_comment', False)
        tree_prefix = options.get('tree_prefix', 'STATE_')

        nexus_keys = {
            'include_translate',
            'include_tree_comment',
            'mode',
            'tree_prefix',
        }
        newick_options = {k: v for k, v in options.items() if k not in nexus_keys}

        if not isinstance(trees, list):
            trees = [trees]

        if isinstance(path_or_stream, str):
            f = open(path_or_stream, mode)
        else:
            f = path_or_stream

        f.write("#NEXUS\n")
        f.write("Begin trees;\n")

        if include_translate:
            newick_options['translator'] = NexusWriter.write_translate_and_map(
                f, trees[0].taxon_names
            )

        tree_counter = 0
        for tree in trees:
            tree_counter += 1
            name = f"{tree_prefix}{tree_counter}"
            newick = tree.newick(**newick_options)
            rooting = 'R' if tree.is_rooted else 'U'
            f.write(f"tree {name} ")
            if include_tree_comment and tree.comment:
                f.write(f"{tree.comment} ")
            f.write(f"= [&{rooting}] {newick}\n")
        f.write("End;\n")

        if isinstance(path_or_stream, str):
            f.close()
