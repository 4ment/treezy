# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

from treezy.bitset import BitSet
from treezy.tree import Tree


def test_create_tree_from_newick():
    newick = "((A:0.1,B:0.2),C:2);"
    taxon_names = ["A", "B", "C"]
    tree = Tree.from_newick(newick, taxon_names)

    assert tree.node_count == 5
    assert tree.leaf_node_count == 3
    assert tree.internal_node_count == 2
    assert tree.root.name is None
    assert tree.is_rooted


def test_compare_taxon_names():
    newick1 = "((A,B),C);"
    newick2 = "((B,C),A);"
    newick3 = "((C,A),B);"

    taxon_names = ["A", "B", "C"]
    taxon_indices = {"A": 0, "B": 1, "C": 2}

    tree0 = Tree.from_newick(newick1)
    tree1 = Tree.from_newick(newick1, taxon_names)
    tree2 = Tree.from_newick(newick2, taxon_names)
    tree3 = Tree.from_newick(newick3, taxon_names)

    assert tree0.taxon_names == taxon_names
    assert tree1.taxon_names == taxon_names

    def check_leaf_ids(tree):
        for node in tree.root.postorder():
            if node.is_leaf:
                assert node.name in taxon_indices
                assert node.id == taxon_indices[node.name]

    check_leaf_ids(tree0)
    check_leaf_ids(tree1)
    check_leaf_ids(tree2)
    check_leaf_ids(tree3)


def test_create_tree_from_newick_with_comments():
    newick = "((A:0.1,B:0.2)[&key=value]:0.1,C:[&a=1]2);"
    taxon_names = ["A", "B", "C"]
    tree = Tree.from_newick(newick, taxon_names)
    tree.root.child_at(0).parse_comment({"key": lambda val: val})
    tree.root.child_at(1).parse_branch_comment({"a": lambda val: int(val)})

    assert tree.node_count == 5
    assert tree.leaf_node_count == 3
    assert tree.internal_node_count == 2
    assert tree.root.name is None
    assert tree.is_rooted
    assert tree.root.child_at(0).annotations["key"]
    assert "a" in tree.root.child_at(1).branch_annotations


def test_reroot_tree():
    newick = "(A:0.1,B:0.2,C:2);"
    taxon_names = ["A", "B", "C"]
    tree = Tree.from_newick(newick, taxon_names)

    assert not tree.is_rooted
    assert tree.make_rooted()
    assert tree.is_rooted
    assert not tree.make_rooted()


def test_make_binary_tree():
    newick = "((A:0.1,B:0.2,C:0.3),D:0.4);"
    taxon_names = ["A", "B", "C", "D"]
    tree = Tree.from_newick(newick, taxon_names)

    assert tree.node_count == 6
    assert tree.leaf_node_count == 4
    assert tree.internal_node_count == 2

    assert tree.make_binary()
    assert tree.node_count == 7


def test_newick_export():
    newick = "((A[&a=b]:0.1,B:[&cat=1]0.2):0.3,C:2);"
    taxon_names = ["A", "B", "C"]
    tree = Tree.from_newick(newick, taxon_names)

    assert tree.newick() == "((A:0.1,B:0.2):0.3,C:2.0);"

    options = {"annotation_keys": ["a"]}
    tree.node_from_id(0).parse_comment()
    assert tree.newick(**options) == "((A[&a=b]:0.1,B:0.2):0.3,C:2.0);"

    options["branch_annotation_keys"] = ["cat"]
    tree.leaf_from_name('B').parse_branch_comment()
    assert tree.newick(**options) == "((A[&a=b]:0.1,B:[&cat=1]0.2):0.3,C:2.0);"

    options = {"annotation_keys": ["key"]}
    tree.root.child_at(0).annotations["key"] = "value"
    assert tree.newick(**options) == "((A:0.1,B:0.2)[&key=value]:0.3,C:2.0);"

    options = {"include_comment": True}
    assert tree.newick(**options) == "((A[&a=b]:0.1,B:0.2):0.3,C:2.0);"

    options = {"decimal_precision": 2}
    assert tree.newick(**options) == "((A:0.10,B:0.20):0.30,C:2.00);"

    newick = "((A,B),C);"
    taxon_names = ["A", "B", "C"]
    tree = Tree.from_newick(newick, taxon_names)

    assert tree.newick() == "((A,B),C);"


def test_compute_descendant_bitset():
    newick = "((A:0.1,B:0.2),C:2);"
    taxon_names = ["A", "B", "C"]
    tree = Tree.from_newick(newick, taxon_names)

    tree.compute_descendant_bitset()
    all_true = ~BitSet(3)

    assert tree.root.descendant_bitset == all_true


def test_post_order_iterator():
    newick = "((A:0.1,B:0.2),C:2);"
    taxon_names = ["A", "B", "C"]
    tree = Tree.from_newick(newick, taxon_names)

    nodes = list(tree.root.postorder())
    assert len(nodes) == 5
    assert nodes[0].name == "A"
    assert nodes[1].name == "B"
    assert nodes[3].name == "C"
    assert nodes[4].is_root


def test_pre_order_iterator():
    newick = "((A:0.1,B:0.2),C:2);"
    taxon_names = ["A", "B", "C"]
    tree = Tree.from_newick(newick, taxon_names)

    nodes = list(tree.root.preorder())
    assert len(nodes) == 5

    assert nodes[0] == tree.root
    assert nodes[2].name == "A"
    assert nodes[3].name == "B"
    assert nodes[4].name == "C"


def test_parse_comment():
    newick = "((A:0.1,B:0.2)[&key=value]:0.1,C:2);"
    taxon_names = ["A", "B", "C"]
    tree = Tree.from_newick(newick, taxon_names)

    converters = {"key": lambda val: val}
    for node in tree.root.postorder():
        node.parse_comment(converters)

    child = tree.root.child_at(0)
    assert "key" not in child.branch_annotations
    assert "key" in child.annotations
    assert child.annotations["key"] == "value"
    del child.annotations["key"]
    assert "key" not in child.annotations
    child.annotations["key"] = "new_value"
    assert child.annotations["key"] == "new_value"

    child.branch_annotations["branch key"] = "value"
    assert "branch key" not in child.annotations


def test_reroot_above():
    newick = "((A:0.1,B:0.2):0.3,C:0.4);"
    taxon_names = ["A", "B", "C"]
    tree = Tree.from_newick(newick, taxon_names)

    tree.reroot_above(tree.root)
    assert tree.newick() == "((A:0.1,B:0.2):0.3,C:0.4);"

    node_to_reroot = tree.root.child_at(0).child_at(0)
    tree.reroot_above(node_to_reroot)
    assert tree.newick() == "(A:0.05,(B:0.2,C:0.7):0.05);"


def test_reroot_above2():
    newick = "(((A:1,B:2):3,C:4):5,D:6);"
    taxon_names = ["A", "B", "C", "D"]
    tree = Tree.from_newick(newick, taxon_names)

    node_to_reroot = tree.root.child_at(0).child_at(0).child_at(0)
    tree.reroot_above(node_to_reroot)
    assert tree.newick() == "(A:0.5,(B:2.0,(C:4.0,D:11.0):3.0):0.5);"

    for node in tree.root.postorder():
        if node != tree.root:
            for child in node.children:
                assert child.parent == node
    assert tree.root.is_root


def test_make_rooted():
    newick = "(A:0.1,B:0.2,C:0.4);"
    taxon_names = ["A", "B", "C"]
    tree = Tree.from_newick(newick, taxon_names)

    node_to_reroot = tree.root.child_at(1)
    assert node_to_reroot.name == "B"
    tree.reroot_above(node_to_reroot)
    assert tree.newick() == "(B:0.1,(A:0.1,C:0.4):0.1);"

    tree = Tree.from_newick(newick, taxon_names)
    tree.make_rooted()
    assert tree.newick() == "(A:0.05,(B:0.2,C:0.4):0.05);"


def test_height_to_distance():
    tree = Tree.from_newick(
        '((A:[&rate=0.1]1,B:[&rate=0.2]2):[&rate=0.1]3,C:[&rate=0.4]4);'
    )
    for node in tree.postorder():
        if not node.is_root:
            node.parse_branch_comment({"rate": lambda rate: float(rate)})
            node.distance *= node.branch_annotations['rate']

    options = {"decimal_precision": 1}
    assert tree.newick(**options) == '((A:0.1,B:0.4):0.3,C:1.6);'
