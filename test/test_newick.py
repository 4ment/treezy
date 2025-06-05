# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

import io

import pytest

from treezy.newick import NewickReader


def test_newickreader_count_trees():
    data = "(A,B);\n(C,D);\n"
    f = io.StringIO(data)
    nf = NewickReader(f)
    assert nf.count_trees() == 2


def test_newickreader_parse():
    taxa_names = ["B", "A"]
    data = "(A,B);\n(B,A);"
    f = io.StringIO(data)
    nf = NewickReader(f, taxa_names)
    trees = nf.parse()
    assert trees[0].taxon_names == taxa_names
    assert trees[1].taxon_names == taxa_names
    assert trees[0].newick() == "(A,B);"
    assert trees[1].newick() == "(B,A);"


def test_newickreader_next_and_has_next():
    taxa_names = []
    data = "(A,B);\n(B,A);"
    f = io.StringIO(data)
    nf = NewickReader(f, taxa_names)
    assert nf.has_next()
    t1 = nf.next()
    assert t1.taxon_names == taxa_names == ["A", "B"]
    assert t1.newick() == "(A,B);"
    assert nf.has_next()
    t2 = nf.next()
    assert t2.newick() == "(B,A);"
    assert t2.taxon_names == taxa_names == ["A", "B"]
    assert not nf.has_next()
    assert nf.next() is None


def test_newickreader_skip_next():
    data = "(A,B);\n(B,A);\n"
    f = io.StringIO(data)
    nf = NewickReader(f)
    assert nf.has_next()
    nf.skip_next()
    assert nf.has_next()
    t = nf.next()
    assert t.newick() == "(B,A);"
    assert not nf.has_next()


def test_newickreader_empty():
    f = io.StringIO("")
    nf = NewickReader(f)
    assert nf.count_trees() == 0
    assert nf.parse() == []
    assert not nf.has_next()
    assert nf.next() is None


def test_newickreader_non_tree_lines():
    data = "notatree\n(A,B);\n#comment\n(B,A);"
    f = io.StringIO(data)
    nf = NewickReader(f)
    assert nf.count_trees() == 2
    trees = nf.parse()
    assert trees[0].newick() == "(A,B);"
    assert trees[1].newick() == "(B,A);"


def test_newickreader_skip_next_on_empty():
    f = io.StringIO("")
    nf = NewickReader(f)
    nf.skip_next()  # Should not raise


def test_newickreader_skip_next_multiple():
    data = "(A,B);\n(A,B);\n(B,A);"
    f = io.StringIO(data)
    nf = NewickReader(f)
    nf.skip_next()
    nf.skip_next()
    assert nf.has_next()
    t = nf.next()
    assert t.newick() == "(B,A);"
    assert not nf.has_next()


def test_newickreader_raise_on_different_taxa():
    data = "(A,B);\n(C,A);"
    f = io.StringIO(data)
    nf = NewickReader(f)
    with pytest.raises(Exception):
        nf.parse()
