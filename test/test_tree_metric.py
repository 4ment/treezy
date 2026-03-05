# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

import pytest

from treezy.tree import Tree
from treezy.tree_metric import RobinsonFouldsMetric


def test_compute_from_sets_identical():
    rf = RobinsonFouldsMetric()
    set1 = {1, 2, 3}
    set2 = {1, 2, 3}
    assert rf.compute_from_splits(set1, set2) == 0.0


def test_compute_from_sets_disjoint():
    rf = RobinsonFouldsMetric()
    set1 = {1, 2}
    set2 = {3, 4}
    assert rf.compute_from_splits(set1, set2) == 4.0


def test_compute_from_sets_partial_overlap():
    rf = RobinsonFouldsMetric()
    set1 = {1, 2, 3}
    set2 = {3, 4}
    # shared = 1, total = 5, so 5 - 2*1 = 3
    assert rf.compute_from_splits(set1, set2) == 3.0


def test_compute_empty_sets():
    rf = RobinsonFouldsMetric()
    set1 = set()
    set2 = set()
    assert rf.compute_from_splits(set1, set2) == 0.0


def test_unrooted():
    tree1 = Tree.from_newick("(A,B,(C,D));")
    tree2 = Tree.from_newick("(C,D,(A,B));", tree1.taxon_names)
    rf = RobinsonFouldsMetric()
    tree1.compute_descendant_bitset()
    tree2.compute_descendant_bitset()
    assert rf.compute(tree1, tree2) == 0.0


def test_rooted_0():
    tree1 = Tree.from_newick("((A,B),(C,D));")
    tree2 = Tree.from_newick("((C,D),(A,B));", tree1.taxon_names)
    tree1.compute_descendant_bitset()
    tree2.compute_descendant_bitset()
    rf = RobinsonFouldsMetric()
    assert rf.compute(tree1, tree2) == 0.0


def test_rooted_2():
    tree1 = Tree.from_newick("((A,B),(C,D));")
    tree2 = Tree.from_newick("(A,(B,(C,D)));", tree1.taxon_names)
    tree1.compute_descendant_bitset()
    tree2.compute_descendant_bitset()
    rf = RobinsonFouldsMetric()
    assert rf.compute(tree1, tree2) == 2.0


def test_compute_weighted_from_dicts_identical():
    rf = RobinsonFouldsMetric()
    splits1 = {frozenset({1, 2}): 0.5, frozenset({3, 4}): 0.3}
    splits2 = {frozenset({1, 2}): 0.5, frozenset({3, 4}): 0.3}
    assert rf.compute_from_splits(splits1, splits2) == 0.0


def test_compute_weighted_from_dicts_different_weights():
    rf = RobinsonFouldsMetric()
    splits1 = {frozenset({1, 2}): 0.5, frozenset({3, 4}): 0.3}
    splits2 = {frozenset({1, 2}): 0.7, frozenset({3, 4}): 0.3}
    # |0.5 - 0.7| + |0.3 - 0.3| = 0.2
    assert rf.compute_from_splits(splits1, splits2) == pytest.approx(0.2)


def test_compute_weighted_from_dicts_disjoint():
    rf = RobinsonFouldsMetric(weighted=True)
    splits1 = {frozenset({1, 2}): 0.5}
    splits2 = {frozenset({3, 4}): 1.5}
    # |0.5 - 0| + |0 - 1.5| = 2.0
    assert rf.compute_from_splits(splits1, splits2) == 2.0


def test_compute_weighted_from_dicts_partial_overlap():
    rf = RobinsonFouldsMetric(weighted=True)
    splits1 = {frozenset({1, 2}): 0.5, frozenset({3, 4}): 0.3}
    splits2 = {frozenset({1, 2}): 0.5, frozenset({5, 6}): 0.2}
    # |0.5 - 0.5| + |0.3 - 0| + |0 - 0.2| = 0.5
    assert rf.compute_from_splits(splits1, splits2) == 0.5


def test_compute_weighted_from_dicts_empty():
    rf = RobinsonFouldsMetric(weighted=True)
    splits1 = {}
    splits2 = {}
    assert rf.compute_from_splits(splits1, splits2) == 0.0


def test_compute_weighted_unrooted_identical():
    tree1 = Tree.from_newick("(A:0.1,B:0.2,(C:0.15,D:0.25):0.35);")
    tree2 = Tree.from_newick("((A:0.1,B:0.2):0.35,C:0.15,D:0.25);", tree1.taxon_names)
    tree1.compute_descendant_bitset()
    tree2.compute_descendant_bitset()
    rf = RobinsonFouldsMetric(weighted=True)
    # # |0.35 - 0.35| = 0.0
    assert rf.compute(tree1, tree2) == 0.0


def test_compute_weighted_unrooted_same_split_different_weights():
    tree1 = Tree.from_newick("(A:0.1,B:0.2,(C:0.15,D:0.25):0.35);")
    tree2 = Tree.from_newick("((A:0.1,B:0.2):0.15,C:0.15,D:0.25);", tree1.taxon_names)
    tree1.compute_descendant_bitset()
    tree2.compute_descendant_bitset()
    rf = RobinsonFouldsMetric(weighted=True)
    # |0.35 - 0.15| = 0.2
    assert rf.compute(tree1, tree2) == pytest.approx(0.2)


def test_compute_weighted_unrooted_distjoint():
    tree1 = Tree.from_newick("(A:0.1,C:0.2,(B:0.15,D:0.25):0.35);")
    tree2 = Tree.from_newick("((A:0.1,B:0.2):0.15,C:0.15,D:0.25);", tree1.taxon_names)
    tree1.compute_descendant_bitset()
    tree2.compute_descendant_bitset()
    rf = RobinsonFouldsMetric(weighted=True)
    # |0.35 - 0.0| + |0 - 0.15| = 0.5
    assert rf.compute(tree1, tree2) == 0.5
