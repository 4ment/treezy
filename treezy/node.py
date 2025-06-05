# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

import weakref
from collections import deque
from io import StringIO
from typing import Any, Dict, Iterator, List, Optional

from treezy.bitset import BitSet


class Node:
    def __init__(self, name: str = None):
        self._name: str = name
        self._id: int = 0
        self._parent: Optional[weakref.ReferenceType['Node']] = None
        self._children: List['Node'] = []
        self._distance: Optional[float] = None
        self._annotations: Dict[str, Any] = {}
        self._branch_annotations: Dict[str, Any] = {}
        self._comment: str = ""
        self._branch_comment: str = ""
        self._descendant_bitset: Optional[BitSet] = None

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, value: str):
        self._name = value

    @property
    def id(self) -> int:
        return self._id

    @id.setter
    def id(self, value: int):
        self._id = value

    @property
    def distance(self) -> Optional[float]:
        return self._distance

    @distance.setter
    def distance(self, value: float):
        self._distance = value

    @property
    def children(self) -> List['Node']:
        return self._children

    def child_count(self) -> int:
        return len(self._children)

    def child_at(self, index: int) -> 'Node':
        return self._children[index]

    def add_child(self, node: 'Node') -> bool:
        if node not in self._children:
            self._children.append(node)
            node.parent = self
            return True
        return False

    def remove_child(self, node: 'Node') -> bool:
        if node in self._children:
            self._children.remove(node)
            node.remove_parent()
            return True
        return False

    def remove_parent(self):
        self._parent = None

    @property
    def parent(self) -> Optional['Node']:
        return self._parent() if self._parent else None

    @parent.setter
    def parent(self, parent: 'Node'):
        self._parent = weakref.ref(parent)

    def siblings(self) -> List['Node']:
        p = self.parent
        if not p:
            return []
        return [n for n in p.children if n is not self]

    @property
    def is_root(self) -> bool:
        return self.parent is None

    @property
    def is_leaf(self) -> bool:
        return len(self._children) == 0

    def contains_annotation(self, key: str) -> bool:
        return key in self._annotations

    def set_annotation(self, key: str, value: Any):
        self._annotations[key] = value

    def annotation(self, key: str) -> Any:
        if key not in self._annotations:
            raise KeyError(f"Key not found: {key}")
        value = self._annotations[key]
        return value

    def contains_branch_annotation(self, key: str) -> bool:
        return key in self._branch_annotations

    def set_branch_annotation(self, key: str, value: Any):
        self._branch_annotations[key] = value

    def branch_annotation(self, key: str) -> Any:
        if key not in self._branch_annotations:
            raise KeyError(f"Key not found: {key}")
        value = self._branch_annotations[key]
        return value

    @property
    def comment(self) -> str:
        return self._comment

    @comment.setter
    def comment(self, comment: str):
        self._comment = comment

    @property
    def branch_comment(self) -> str:
        return self._branch_comment

    @branch_comment.setter
    def branch_comment(self, comment: str):
        self._branch_comment = comment

    def remove_annotation(self, key: str):
        if key in self._annotations:
            del self._annotations[key]

    def clear_annotations(self):
        self._annotations.clear()

    def annotation_keys(self) -> List[str]:
        return list(self._annotations.keys())

    def make_binary(self) -> bool:
        made_binary = False
        while self.child_count() > 2:
            child0 = self._children[0]
            child1 = self._children[1]
            self.remove_child(child0)
            self.remove_child(child1)
            new_node = Node()
            new_node.add_child(child0)
            new_node.add_child(child1)
            new_node.distance = 0
            self._children.insert(0, new_node)
            new_node.parent = self
            made_binary = True
        return made_binary

    def is_binary(self) -> bool:
        return len(self._children) == 2

    @property
    def descendant_bitset(self) -> 'BitSet':
        if not hasattr(self, '_descendant_bitset'):
            raise AttributeError(
                "BitSet not computed. Call compute_descendant_bitset(size) first."
            )
        return self._descendant_bitset

    def compute_descendant_bitset(self, size) -> None:
        if self.is_leaf:
            self._descendant_bitset = BitSet(size)
            self._descendant_bitset[self._id] = True
        else:
            self._descendant_bitset = BitSet(size)
            for child in self.children:
                self._descendant_bitset |= child._descendant_bitset

    def newick(self, options: Optional[Dict[str, Any]] = None) -> str:
        """Generate a Newick string representation of the node.

        This method constructs a Newick format string for the node and its children,
        including branch lengths and comments if specified in the options.
        The options dictionary can include:
            - include_branch_lengths (bool): Whether to include branch lengths
            in the output.
            - decimal_precision (int): Number of decimal places for branch lengths.
            - include_internal_node_name (bool): Whether to include internal node names.
            - translator (Dict[str, str]): A mapping for translating node names.
        If `translator` is provided, it will be used to translate node names
        in the output.

        Args:
            options (Optional[Dict[str, Any]], optional): A dictionary of options for
            formatting the Newick string. If None, defaults will be used.

        Returns:
            str: A Newick formatted string representing the node and its children.
        """
        options = options or {}
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
            children_strs = [child.newick(options) for child in self.children]
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

    def parse_comment(self, converters: Dict[str, Any] = None) -> Dict[str, Any]:
        if not self._comment:
            return {}
        annotations = parse_comment(self._comment, converters)
        if annotations:
            self._annotations.update(annotations)
        return annotations

    def parse_branch_comment(self, converters: Dict[str, Any] = None) -> None:
        if not self._branch_comment:
            return {}
        annotations = parse_comment(self._branch_comment, converters)
        if annotations:
            self._branch_annotations.update(annotations)

    def _make_comment_for_newick(self, options: Dict[str, Any]) -> tuple[str, str]:
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
            self._comment, self._annotations, "include_comment", "annotation_keys"
        )
        branch_comment = build_comment(
            self._branch_comment,
            self._branch_annotations,
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
        stack = [self]
        out = deque()  # acts like a reverse postorder collector

        while stack:
            node = stack.pop()
            out.appendleft(node)  # reverse insertion
            stack.extend(node.children)  # children left-to-right

        return iter(out)

    def preorder(self) -> Iterator['Node']:
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))

    def levelorder(self) -> Iterator['Node']:
        queue = deque([self])
        while queue:
            node = queue.popleft()
            yield node
            queue.extend(node.children)

    def __repr__(self) -> str:
        return f"Node(name={self._name}, id={self._id}, distance={self._distance})"


def parse_comment(comment: str, converters: Dict[str, Any] = None):
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
