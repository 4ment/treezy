# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

import copy
import weakref
from collections import deque
from io import StringIO
from typing import Any, Dict, Iterator, List, Optional

from treezy.bitset import BitSet


class Node:
    """Represents a node in a tree-like structure.

    Each node may have a name, an ID, a distance to its parent, child nodes, and
    a parent node. It can contain annotations and comments associated with the node
    or its branch leading to the parent.

    In practice:

    - Leaf nodes (i.e., nodes without children) must have a name.
    - Internal nodes are not required to have a name.
    - The `id` can be arbitrary for standalone nodes, but within a `Tree` object it
      is unique,ranging from 0 to the total number of nodes.
    - Leaf node IDs are guaranteed to be less than the number of taxa.
    - The `distance` field may be omitted if only the topology is of interest.

    Attributes
    ----------
    name
        The name of the node.
    id
        Unique identifier for the node.
    distance
        distance to the parent node.
    children
        Child nodes of this node.
    annotations
        Node-level annotations.
    branch_annotations
        Annotations on the branch to the parent.
    comment
        comment associated with the node.
    branch_comment
        comment associated with the branch to parent.
    """

    name: Optional[str]
    id: int
    _parent: Optional[weakref.ReferenceType['Node']]
    children: List['Node']
    distance: Optional[float]
    annotations: Dict[str, Any]
    branch_annotations: Dict[str, Any]
    comment: Optional[str]
    branch_comment: Optional[str]
    _descendant_bitset: Optional[BitSet]

    def __init__(self, name: Optional[str] = None):
        """Initializes a Node with an optional name.

        Parameters
        ----------
        name
            the name of the node. Defaults to None.
        """
        self.name = name
        self.id = -1  # a sensible id will be set later by a Tree object
        self._parent = None
        self.children = []
        self.distance = None
        self.annotations = {}
        self.branch_annotations = {}
        self.comment = None
        self.branch_comment = None
        self._descendant_bitset = None

    def __deepcopy__(self, memo):
        copied = Node(copy.deepcopy(self.name, memo))
        memo[id(self)] = copied
        copied.id = self.id
        copied.distance = copy.deepcopy(self.distance, memo)
        copied.annotations = copy.deepcopy(self.annotations, memo)
        copied.branch_annotations = copy.deepcopy(self.branch_annotations, memo)
        copied.comment = self.comment
        copied.branch_comment = self.branch_comment
        copied._descendant_bitset = copy.deepcopy(self._descendant_bitset, memo)

        copied.children = [copy.deepcopy(child, memo) for child in self.children]
        for child in copied.children:
            child._parent = copied

        return copied

    def __getitem__(self, index: int) -> 'Node':
        return self.children[index]

    def __iter__(self):
        return iter(self.children)

    def __contains__(self, node: 'Node') -> bool:
        return node in self.children

    def child_at(self, index: int) -> 'Node':
        """Get the child node at the specified index.

        Parameters
        ----------
        index
            Index of the child node.

        Returns
        -------
        Node
            The child node at the specified index.
        """
        return self.children[index]

    def add_child(self, node: 'Node') -> bool:
        """Add a child node to this node.

        Parameters
        ----------
        node
            The child node to add.

        Returns
        -------
        bool
            True if the child was added, False if it was already present.
        """
        if node not in self.children:
            self.children.append(node)
            node.parent = self
            return True
        return False

    def remove_child(self, node: 'Node') -> bool:
        """Remove a child node from this node.

        Parameters
        ----------
        node
            The child node to remove.

        Returns
        -------
        bool
            True if the child was removed, False if it was not found.
        """
        if node in self.children:
            self.children.remove(node)
            node.parent = None
            return True
        return False

    @property
    def parent(self) -> Optional['Node']:
        """Get or set the parent node of this node."""
        return self._parent() if self._parent else None

    @parent.setter
    def parent(self, parent: 'Node'):
        self._parent = weakref.ref(parent) if parent else None

    def siblings(self) -> List['Node']:
        """Get a list of sibling nodes."""
        p = self.parent
        if not p:
            return []
        return [n for n in p.children if n is not self]

    @property
    def is_root(self) -> bool:
        """Check if this node is the root of the tree.

        A root node is defined as a node that has no parent.

        Returns
        -------
        bool
            True if this node is the root, False otherwise.
        """
        return self.parent is None

    @property
    def is_leaf(self) -> bool:
        """Check if this node is a leaf node."""
        return len(self.children) == 0

    def collapse(self) -> None:
        """Collapse this node into its parent.

        This method removes this node from its parent's children and adds all
        of its children to the parent. After this operation, this node will no
        longer be part of the tree, and its children will be directly under
        its parent.
        """
        parent = self.parent
        parent.remove_child(self)
        for child in self.children:
            parent.add_child(child)
        self.children.clear()

    def make_binary(self) -> bool:
        """Convert this node into a binary node if it has more than two children.
        If the node has more than two children, it will create new internal nodes
        to ensure that each internal node has at most two children. The new nodes
        will have a distance of 0, and the original children will be moved under
        these new nodes.

        Returns
        -------
        bool
            True if the node was made binary, False if it was already binary.
        """
        made_binary = False
        while len(self.children) > 2:
            child0 = self.children[0]
            child1 = self.children[1]
            self.remove_child(child0)
            self.remove_child(child1)
            new_node = Node()
            new_node.add_child(child0)
            new_node.add_child(child1)
            new_node.distance = 0
            self.children.insert(0, new_node)
            new_node.parent = self
            made_binary = True
        return made_binary

    def is_binary(self) -> bool:
        """Check if this node is binary, meaning it has at most two children.

        Returns
        -------
        bool
            True if the node is binary (has 2 children), False otherwise.
        """
        return len(self.children) == 2

    @property
    def descendant_bitset(self) -> 'BitSet':
        """Get the BitSet representing all descendants of this node.

        Returns
        -------
        BitSet
            A bitset where each bit represents whether a node is a descendant
            of this node.

        Raises
        ------
        AttributeError
            If the BitSet has not been computed yet.
        """
        if self._descendant_bitset is None:
            raise AttributeError(
                "BitSet not computed. Call compute_descendant_bitset(size) first."
            )
        return self._descendant_bitset

    def compute_descendant_bitset(self, size: int) -> None:
        """Compute the BitSet of all descendants of this node.
        This method populates the `_descendant_bitset` attribute with a BitSet
        that represents all nodes that are descendants of this node.
        If the node is a leaf, it will set the bit corresponding to its own ID.
        If the node has children, it will recursively compute the BitSet for each child
        and combine them into a single BitSet for this node.

        Parameters
        ----------
        size
            The size of the BitSet, which should be the number of taxa in the tree.
        """
        if self.is_leaf:
            self._descendant_bitset = BitSet(size)
            self._descendant_bitset[self.id] = True
        else:
            self._descendant_bitset = BitSet(size)
            for child in self.children:
                self._descendant_bitset |= child.descendant_bitset

    def newick(self, **options) -> str:
        """Generate a Newick string representation of the node.

        This method constructs a Newick format string for the node and its children,
        including branch lengths and comments if specified in the options.
        The options dictionary can include:

            - `include_branch_lengths` (bool): Whether to include branch lengths
              in the output.
            - `decimal_precision` (int): Number of decimal places for branch lengths.
            - `include_internal_node_name` (bool): Whether to include internal
              node names.
            - `translator` (Dict[str, str]): A mapping for translating node names.

        If `translator` is provided, it will be used to translate node names
        in the output.

        Parameters
        ----------
        options
            A dictionary of options for formatting the Newick string.
            If None, defaults will be used.

        Returns
        -------
        str
            A Newick formatted string representing the node and its children.
        """
        include_branch_lengths = options.get("include_branch_lengths", True)
        decimal_precision = options.get("decimal_precision", -1)
        include_internal_node_name = options.get("include_internal_node_name", False)
        translator = options.get("translator", None)

        def format_distance(dist):
            if decimal_precision > 0:
                return f"{dist:.{decimal_precision}f}"
            else:
                return f"{dist}"

        buf = StringIO()
        if self.is_leaf:
            if translator and self.name in translator:
                buf.write(translator[self.name])
            else:
                buf.write(self.name)
            comment, branch_comment = self._make_comment_for_newick(options)
            if include_branch_lengths and self.distance is not None:
                buf.write(comment)
                buf.write(":")
                buf.write(branch_comment)
                buf.write(format_distance(self.distance))
            else:
                buf.write(comment)
        else:
            buf.write("(")
            children_strs = [child.newick(**options) for child in self.children]
            buf.write(",".join(children_strs))
            buf.write(")")
            comment, branch_comment = self._make_comment_for_newick(options)
            if include_internal_node_name and self.name:
                buf.write(self.name)
            if include_branch_lengths and self.distance is not None:
                buf.write(comment)
                buf.write(":")
                buf.write(branch_comment)
                buf.write(format_distance(self.distance))
            else:
                buf.write(comment)
        return buf.getvalue()

    def parse_comment(self, converters: Dict[str, Any] = None) -> None:
        """Parse the comment associated with this node.

        This method extracts annotations from the comment string using the provided
        converters dictionary, which maps annotation keys to conversion functions.

        Example
        -------
        If the comment is ``[&rate=0.01,name=myname]``, and `converters` is::

            {"rate": lambda x: float(x)}

        then this method will return::

            {"rate": 0.01, "name": "myname"}

        Parameters
        ----------
        converters
            A dictionary mapping annotation keys to conversion functions.
            Defaults to None.
        """
        if self.comment is not None:
            annotations = parse_comment(self.comment, converters)
            if annotations:
                self.annotations.update(annotations)

    def parse_branch_comment(self, converters: Dict[str, Any] = None) -> None:
        """Parse the branch comment associated with this node.

        This method extracts annotations from the comment string using the provided
        converters dictionary, which maps annotation keys to conversion functions.

        Example
        -------
        If the comment is ``[&rate=0.01,name=myname]``, and `converters` is::

            {"rate": lambda x: float(x)}

        then this method will return::

            {"rate": 0.01, "name": "myname"}

        Parameters
        ----------
        converters
            A dictionary mapping annotation keys to conversion functions.
            Defaults to None.
        """
        if self.branch_comment is not None:
            annotations = parse_comment(self.branch_comment, converters)
            if annotations:
                self.branch_annotations.update(annotations)

    def _make_comment_for_newick(self, options: Dict[str, Any]) -> tuple[str, str]:
        """Generate comments for Newick output based on options.

        This method builds comments for the node and its branch based on the
        provided options.

        Parameters
        ----------
        options
            A dictionary of options that may include:
                - include_comment (bool): Whether to include the node comment.
                - include_branch_comment (bool): Whether to include the branch comment.
                - annotation_keys (list): List of keys to include in the node comment.
                - branch_annotation_keys (list): List of keys to include in the branch
                  comment.

        Returns
        -------
        tuple[str, str]
            a tuple containing the node comment and branch comment.
        """

        def build_comment(raw_comment, annotations, include_flag, keys_option):
            if options.get(include_flag, False) and raw_comment is not None:
                return raw_comment
            elif keys_option in options:
                comment = "[&"
                for key in options[keys_option]:
                    if key in annotations:
                        val = annotations[key]
                        comment += f"{key}={val},"
                if comment.endswith(","):
                    comment = comment[:-1]
                comment += "]"
                if comment == "[&]":
                    comment = ""
                return comment
            return ""

        comment = build_comment(
            self.comment, self.annotations, "include_comment", "annotation_keys"
        )
        branch_comment = build_comment(
            self.branch_comment,
            self.branch_annotations,
            "include_branch_comment",
            "branch_annotation_keys",
        )
        return comment, branch_comment

    # def postorder(self) -> Iterator['Node']:
    #     stack = deque([self])
    #     while stack:
    #         node = stack.pop()
    #         yield node
    #         stack.extend(node.children)

    def postorder(self) -> Iterator['Node']:
        """Generate nodes in postorder traversal.

        This method yields nodes in postorder, meaning it visits all children
        before the parent node.

        Returns
        -------
        Iterator[Node]
            An iterator that yields nodes in postorder.
        """
        stack = [self]
        out = deque()  # acts like a reverse postorder collector

        while stack:
            node = stack.pop()
            out.appendleft(node)  # reverse insertion
            stack.extend(node.children)  # children left-to-right

        return iter(out)

    def preorder(self) -> Iterator['Node']:
        """Generate nodes in preorder traversal.

        This method yields nodes in preorder, meaning it visits the parent node
        before its children.

        Returns
        -------
        Iterator[Node]
            An iterator that yields nodes in preorder.
        """
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))

    def levelorder(self) -> Iterator['Node']:
        """Generate nodes in level order traversal.

        This method yields nodes in level order, meaning it visits all nodes at
        the current level before moving to the next level.

        Returns
        -------
        Iterator[Node]
            An iterator that yields nodes in level order.
        """
        queue = deque([self])
        while queue:
            node = queue.popleft()
            yield node
            queue.extend(node.children)

    def __repr__(self) -> str:
        return f"Node(name={self.name}, id={self.id}, distance={self.distance})"


def parse_comment(comment: str, converters: Dict[str, Any] = None):
    """Parse a comment string into a dictionary of annotations.

    This function extracts key-value pairs from a comment string formatted as
    ``[&mean=0.2,hpd={0.1,0.6}]``. It supports nested structures and can
    handle commas inside brackets. If a key is found in the `converters`
    dictionary, the corresponding value will be converted using the provided
    converter function. If no converters are provided, the values will be
    stored as strings.
    The comment string should start with ``&`` and end with ``]``, with key-value
    pairs separated by commas.

    This type of comment is used by BEAST package.

    Parameters
    ----------
    comment
        The comment string to parse, formatted as `[&key1=value1,key2=value2,...]`.
    converters
        A dictionary mapping annotation keys to conversion functions. If a key
        is found in this dictionary, its value will be converted using the
        corresponding function. If None, values will be stored as strings.
    """
    annotations: Dict[str, Any] = {}

    start = comment.find("&") + 1
    end = comment.rfind("]")
    if start == 0 or end == -1 or end <= start:
        return
    content = comment[start:end]

    def split_outside_brackets(s: str) -> list[str]:
        result = []
        buffer = []
        stack = []

        for char in s:
            if char in '[{':
                stack.append(char)
            elif char in ']}':
                if stack and (
                    (char == ']' and stack[-1] == '[')
                    or (char == '}' and stack[-1] == '{')
                ):
                    stack.pop()
            if char == ',' and not stack:
                result.append(''.join(buffer).strip())
                buffer = []
            else:
                buffer.append(char)

        if buffer:
            result.append(''.join(buffer).strip())
        return result

    tokens = split_outside_brackets(content)
    for token in tokens:
        eq_pos = token.find('=')
        if eq_pos != -1:
            key = token[:eq_pos].strip()
            value = token[eq_pos + 1 :].strip()
            if converters and key in converters:
                annotations[key] = converters[key](value)
            else:
                annotations[key] = value

    return annotations
