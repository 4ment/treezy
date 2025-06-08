# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

import copy
from random import choice
from typing import Any, Dict, Iterator, List, Optional

# from treezy.newick import tokenize_newick
from treezy.node import Node


class Tree:
    """Represents a phylogenetic tree.

    Provides methods to create, manipulate, and traverse a tree structure.
    Supports operations such as rerooting, binarization, and Newick export.
    Each tree has a root node and a list of taxon names. The order of these
    names corresponds to the identifiers of the leaf nodes.

    Attributes
    ----------
    name
        The name of the tree.
    comment
        A comment associated with the tree. Defaults to None.
    annotations
        Annotations for the tree.
    """

    name: str
    _root: Node
    _taxon_names: Optional[List[str]]
    _nodes: List[Node]
    _leaf_count: int
    _internal_count: int
    _node_count: int
    comment: Optional[str]
    annotations: Dict[str, Any]

    def __init__(self, root: Node, taxon_names: Optional[List[str]] = None):
        """Initialize a Tree object with a root node and optional list taxon names.

        If taxon names are not provided, the tree will extract them from the leaf nodes
        using a postorder traversal.

        Parameters
        ----------
        root
            The root node of the tree.
        taxon_names
            A list of taxon names. If not provided, the tree will extract taxon names
            from the leaf nodes using a postorder traversal. If provided, it will
            validate that the names match those found in the tree structure.
        """
        self.name = None
        self._root = root
        self._taxon_names = taxon_names or []
        self._nodes = []
        self._leaf_count = 0
        self._internal_count = 0
        self._node_count = 0
        self.comment = None
        self.annotations = {}

        if taxon_names is None or len(taxon_names) == 0:
            for node in self._root.postorder():
                if node.is_leaf:
                    self._taxon_names.append(node.name)
                    node.id = self._leaf_count
                    self._leaf_count += 1
            self._nodes = [None] * self._leaf_count

            for node in self._root.postorder():
                if not node.is_leaf:
                    node.id = self._leaf_count + self._internal_count
                    self._internal_count += 1
                    self._nodes.append(node)
                else:
                    self._nodes[node.id] = node

            self._taxon_names = [node.name for node in self._nodes if node.is_leaf]
            self._node_count = self._leaf_count + self._internal_count
        else:
            self.update_ids()

    def __deepcopy__(self, memo):
        copied_root = copy.deepcopy(self._root, memo)
        copied_tree = Tree(copied_root, copy.deepcopy(self._taxon_names, memo))
        memo[id(self)] = copied_tree
        copied_tree.comment = copy.deepcopy(self.comment, memo)
        copied_tree.annotations = copy.deepcopy(self.annotations, memo)

        return copied_tree

    @property
    def taxon_names(self) -> List[str]:
        """Get or set the list of taxon names for the tree."""
        return self._taxon_names

    @taxon_names.setter
    def taxon_names(self, names: List[str]):
        self._taxon_names = names
        self.update_ids()

    @property
    def node_count(self) -> int:
        """Get the total number of nodes in the tree."""
        return self._node_count

    @property
    def internal_node_count(self) -> int:
        """Get the number of internal nodes in the tree."""
        return self._internal_count

    @property
    def leaf_node_count(self) -> int:
        """Get the number of leaf nodes in the tree."""
        return self._leaf_count

    @property
    def root(self) -> Node:
        """Get the root node of the tree."""
        return self._root

    @property
    def nodes(self) -> List[Node]:
        """Get a list of all nodes in the tree."""
        return self._nodes

    def node_from_id(self, node_id: int) -> Node:
        """Get a node by its ID.

        Parameters
        ----------
        node_id
            The ID of the node to retrieve.

        Returns
        -------
        Node
            The node with the specified ID.
        """
        return self._nodes[node_id]

    def leaf_from_name(self, name: str) -> Node:
        """Get a leaf node by its name.

        Parameters
        ----------
        name
            The name of the leaf node to retrieve.

        Returns
        -------
        Node
            The leaf node with the specified name.

        Raises
        ------
        ValueError
            If the name does not exist in taxon_names.
        """
        return self._nodes[self._taxon_names.index(name)]

    def compute_descendant_bitset(self) -> None:
        """Compute in postorder the descendant bitset for each node in the tree."""
        for node in self._root.postorder():
            node.compute_descendant_bitset(self._leaf_count)

    @property
    def is_rooted(self) -> bool:
        """Check if the tree is rooted."""
        return len(self._root.children) == 2

    def make_rooted(self) -> bool:
        """Make the tree rooted if it is not already.

        Returns
        -------
        bool
        True if the tree was not rooted and has been made rooted, False otherwise.
        """
        degree = len(self._root.children)
        if degree > 2:
            self.reroot_above(self._root.children[0])
        return degree > 2

    def make_binary(self) -> bool:
        """Convert the tree to a binary tree if it is not already.

        Returns
        -------
        bool:
        True if the tree was not binary and has been made binary, False otherwise.
        """
        made_binary = False
        for node in reversed(self._nodes):
            if node.is_leaf:
                continue
            if len(node.children) > 2:
                made_binary |= node.make_binary()
        if made_binary:
            self.update_ids()
        return made_binary

    def reroot_above(self, node: Node):
        """Reroot the tree above the specified node.

        This method reroots the tree by creating a new root node between the specified
        node and its parent. It adjusts the distances of the nodes and ensures that the
        tree remains valid. If the specified node is already the root, it does nothing.

        Parameters
        ----------
        node
            The node above which to reroot the tree.
        """
        # If node is already the root, do nothing
        if node.is_root:
            return

        parent = node.parent

        # If node is just below the root
        if parent.is_root:
            midpoint = node.distance / 2
            node.distance -= midpoint

            siblings = node.siblings()

            if len(parent.children) > 2:
                # Make parent binary by grouping siblings under a new node
                new_node = Node()
                for sibling in siblings:
                    parent.remove_child(sibling)
                    new_node.add_child(sibling)
                parent.add_child(new_node)
                new_node.distance = midpoint
                self.update_ids()
            else:
                siblings[0].distance += midpoint
        else:
            new_root = Node()
            grandparent = parent.parent

            new_root.add_child(node)
            new_root.add_child(parent)

            branch_length = parent.distance
            midpoint = node.distance / 2

            node.distance = midpoint
            parent.distance = midpoint

            n = parent
            n_parent = grandparent

            n_parent.remove_child(parent)
            parent.remove_child(node)

            node.parent = new_root
            parent.parent = new_root

            while not n_parent.is_root:
                temp = n_parent.parent
                n.add_child(n_parent)

                bl = n_parent.distance
                n_parent.distance = branch_length
                branch_length = bl

                n_parent.parent = n

                n = n_parent
                n_parent = temp

                if n_parent:
                    temp2 = n.parent
                    n_parent.remove_child(n)
                    n.parent = temp2

            # n_parent is the old root
            if n_parent and n_parent.children:
                unaffected_sibling = n_parent.children[0]
                unaffected_sibling.distance += branch_length
                n.add_child(unaffected_sibling)

            self._root = new_root
            self.update_ids()

    def update_ids(self):
        """Update the IDs of the nodes in the tree based on the taxon names.

        This method assigns IDs to the nodes in the tree based on their position in
        the taxon names list. Leaf nodes are assigned IDs corresponding to their index
        in the taxon names list, while internal nodes are assigned IDs starting from
        the number of leaf nodes. The internal node IDs are assigned in a postorder
        traversal of the tree.
        The method also initializes the `_nodes` list, which holds references to all
        nodes in the tree, allowing for efficient access by ID or we need to iterate
        through all nodes in the tree without a particular order.
        """
        self._internal_count = 0
        self._leaf_count = len(self._taxon_names)
        self._nodes = [None] * self._leaf_count

        for node in self._root.postorder():
            if not node.is_leaf:
                node.id = self._leaf_count + self._internal_count
                self._nodes.append(node)
                self._internal_count += 1
            else:
                node.id = self._taxon_names.index(node.name)
                self._nodes[node.id] = node
        self._node_count = self._leaf_count + self._internal_count

    def newick(self, **options) -> str:
        """Export the tree in Newick format.

        By default, the Newick string will not include comments, branch
        comments, internal node names. You can customize the output by providing
        options in the `options` dictionary.
        For details on the supported options, see the documentation for
        :meth:`Node.newick <treezy.node.Node.newick>`.

        Parameters
        ----------
        options
            A dictionary of options for formatting the Newick string.
            See :meth:`Node.newick <treezy.node.Node.newick>` for supported keys.

        Returns
        -------
        str
            A Newick formatted string representing the tree structure.
        """
        return self._root.newick(**options) + ";"

    @classmethod
    def from_newick(
        cls, newick: str, taxon_names: Optional[List[str]] = None, **options
    ) -> 'Tree':
        """Parse a Newick string to create a Tree object.

        This method parses a Newick formatted string and constructs a tree structure.
        It handles comments, branch lengths, and taxon names. If taxon names are
        provided, it checks that they match the names found in the Newick string.
        If not provided, it extracts taxon names from the Newick string.
        The Newick format is expected to be well-formed, and the method will raise
        a ValueError if the provided taxon names do not match those found in the
        Newick string.
        The Newick string should be in the format: ``(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);``
        where `A`, `B`, `C`, and `D` are taxon names, and the numbers represent
        branch lengths.

        The method supports **node comments** in the Newick string, which can be
        included in square brackets after the taxon names and internal node names
        or closing brackets, like so:

        ``(A[&key=1]:0.1,B[&key=2]:0.2,(C:0.3,D:0.4)[&key=3]:0.5);``

        The Newick string can also include **branch comments**, which are enclosed in
        square brackets after the colon following the taxon name or branch length,
        like so:
        ``(A:0.1[&key=1],B:0.2[&key=2],(C:0.3,D:0.4):[&key3]0.5);``

        Branch and node comments can be used together, like so:
        ``(A[&taxon=A]:0.1[&branch=1],B:0.2,(C:0.3,D:0.4):0.5);``
        which should be interpreted as:
        node A has a node comment `taxon=A` and a branch leading to the root with
        comment `branch=1` with length 0.1.

        The method will extract these comments and store them as strings in the
        corresponding Node objects. By default, it won't parse these comments,
        but you can parse them into structured annotations by calling the
        :meth:`parse_comment<treezy.node.Node.parse_comment>` and
        :meth:`parse_branch_comment<treezy.node.Node.parse_branch_comment>`
        methods on the Node objects after creating the tree.

        If the Newick string is malformed or contains unexpected characters, it will
        raise a ValueError.

        Parameters
        ----------
        newick
            A Newick formatted string representing a tree structure.
        taxon_names
            A list of taxon names to validate against the names found in the Newick
            string. If provided, it will check that the names match those in the Newick
            string. If not provided, it will extract taxon names from the Newick string.
        **kwargs
            Additional keyword arguments. Currently supports:
                - `strip_quotes` (bool): If True, will strip quotes from taxon names
                  in the Newick string. Defaults to False.

        Raises
        ------
        ValueError
            If the provided taxon names do not match those found in the
            Newick string, or if the Newick string is malformed.

        Returns
        -------
        Tree
            A Tree object representing the parsed Newick structure.
        """
        newick = newick.strip()
        stack = []
        current_taxon_names = []
        taxon_counter = 0
        just_closed = False
        i = 0
        n = len(newick)
        strip_quotes = options.get('strip_quotes', False)

        while i < n:
            c = newick[i]
            if c == '[':
                # node comment
                start = i
                while i < n and newick[i] != ']':
                    i += 1
                comment = newick[start : i + 1]
                stack[-1].comment = comment
            elif c == ':':
                i += 1
                start = i
                # branch comment
                if newick[i] == '[':
                    while newick[i] != ']':
                        i += 1
                    i += 1
                    stack[-1].branch_comment = newick[start:i]
                    start = i

                while i < n and (newick[i].isdigit() or newick[i] in '.eE-+'):
                    i += 1
                branch_length_str = newick[start:i]
                i -= 1
                try:
                    branch_length = float(branch_length_str)
                except ValueError:
                    branch_length = 0.0
                stack[-1].distance = branch_length
            elif c not in '(),;':
                start = i
                while i < n and newick[i] not in ':[,);':
                    i += 1
                identifier = newick[start:i].strip()
                i -= 1

                if just_closed:
                    stack[-1].name = identifier
                else:
                    if strip_quotes:
                        identifier = identifier.strip('"\'')
                    node = Node(identifier)
                    node.id = taxon_counter
                    taxon_counter += 1
                    current_taxon_names.append(identifier)
                    stack[-1].add_child(node)
                    stack.append(node)
            elif c == '(':
                just_closed = False
                node = Node()
                if len(stack) != 0:
                    stack[-1].add_child(node)
                stack.append(node)
            elif c in '),':
                stack.pop()
                just_closed = c == ')'
            # else:
            #     print(f"Unexpected character: {c} at position {i}")

            i += 1

        if taxon_names is None:
            taxon_names = current_taxon_names
        elif len(taxon_names) == 0:
            taxon_names.extend(current_taxon_names)
        elif set(taxon_names) != set(current_taxon_names):
            missing = set(current_taxon_names) - set(taxon_names)
            extra = set(taxon_names) - set(current_taxon_names)
            if missing:
                print(f"Missing taxon names: {missing}")
            if extra:
                print(f"Extra taxon names: {extra}")
            raise ValueError(
                "Provided taxon names do not match those in the Newick string."
            )

        root = stack[0] if stack else None
        return cls(root, taxon_names)

    @classmethod
    def random(cls, taxon_names: List[str]) -> 'Tree':
        """Generate a random binary tree from a list of taxon names.
        The tree is generated by randomly combining nodes until only one node remains.


        Parameters
        ----------
        taxon_names
            List of taxon names to create the tree from.

        Returns
        -------
        Tree
            A Tree object with a randomly generated structure based on the provided
            taxon names.
        """
        nodes = [Node(name) for name in taxon_names]
        while len(nodes) > 1:
            index1 = choice(range(len(nodes)))
            index2 = choice(range(len(nodes)))
            while index1 == index2:
                index2 = choice(range(len(nodes)))
            new_node = Node()
            new_node.add_child(nodes[index1])
            new_node.add_child(nodes[index2])
            nodes.append(new_node)
            nodes.pop(max(index1, index2))
        return cls(nodes[0], taxon_names)

    def postorder(self) -> Iterator['Node']:
        """Return an iterator for postorder traversal of the tree."""
        return self._root.postorder()

    def preorder(self) -> Iterator['Node']:
        """Return an iterator for preorder traversal of the tree."""
        return self._root.preorder()

    def levelorder(self) -> Iterator['Node']:
        """Return an iterator for levelorder traversal of the tree."""
        return self._root.levelorder()

    def __str__(self):
        return self.newick()

    def __repr__(self):
        return f"Tree(root={self._root}, taxon_names={self._taxon_names})"

    def __len__(self):
        """Return the number of nodes in the tree."""
        return self.node_count

    def __getitem__(self, index: int) -> Node:
        """Get a node by its index."""
        if index < 0 or index >= self.node_count:
            raise IndexError("Index out of range")
        return self.node_from_id(index)

    def __iter__(self) -> Iterator[Node]:
        """Return an iterator for the tree's nodes."""
        return iter(self._nodes)
