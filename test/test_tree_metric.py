# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

from treezy.tree_metric import RobinsonFouldsMetric


class DummyNode:
    def __init__(self, descendant_bitset, is_leaf=False, is_root=False):
        self.descendant_bitset = descendant_bitset
        self.is_leaf = is_leaf
        self.is_root = is_root


class DummyTree:
    def __init__(self, nodes):
        self.nodes = nodes


def test_compute_from_sets_identical():
    rf = RobinsonFouldsMetric()
    set1 = {1, 2, 3}
    set2 = {1, 2, 3}
    assert rf.compute_from_sets(set1, set2) == 0.0


def test_compute_from_sets_disjoint():
    rf = RobinsonFouldsMetric()
    set1 = {1, 2}
    set2 = {3, 4}
    assert rf.compute_from_sets(set1, set2) == 4.0


def test_compute_from_sets_partial_overlap():
    rf = RobinsonFouldsMetric()
    set1 = {1, 2, 3}
    set2 = {3, 4}
    # shared = 1, total = 5, so 5 - 2*1 = 3
    assert rf.compute_from_sets(set1, set2) == 3.0


def test_compute_empty_sets():
    rf = RobinsonFouldsMetric()
    set1 = set()
    set2 = set()
    assert rf.compute_from_sets(set1, set2) == 0.0


def test_compute_with_trees():
    # Only internal nodes (not leaf, not root) are considered
    nodes1 = [
        DummyNode(1, is_leaf=True),
        DummyNode(2, is_root=True),
        DummyNode(3),
        DummyNode(4),
    ]
    nodes2 = [
        DummyNode(1, is_leaf=True),
        DummyNode(2, is_root=True),
        DummyNode(4),
        DummyNode(5),
    ]
    tree1 = DummyTree(nodes1)
    tree2 = DummyTree(nodes2)
    rf = RobinsonFouldsMetric()
    # Internal nodes: tree1 = {3,4}, tree2 = {4,5}
    # shared = {4}, total = 4, so 4 - 2*1 = 2
    assert rf.compute(tree1, tree2) == 2.0
