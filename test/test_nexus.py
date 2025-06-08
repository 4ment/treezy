# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

import io
import textwrap

import pytest

from treezy.nexus import NexusReader, NexusWriter
from treezy.tree import Tree


def make_nexus(trees, translate=None, taxa=None, comments=None):
    lines = ["#NEXUS", "Begin trees;"]
    if translate:
        lines.append("  Translate")
        translate_map = {}
        for i, (sh, name) in enumerate(translate):
            end = "," if i < len(translate) - 1 else ";"
            lines.append(f"    {sh} {name}{end}")
            translate_map[name] = sh
    for i, t in enumerate(trees):
        if comments:
            lines.append(f"  tree tree{i+1} {comments[i]}= {t}")
        else:
            lines.append(f"  tree tree{i+1} = {t}")
    lines.append("End;")
    if taxa:
        lines.insert(1, "Begin taxa;")
        lines.insert(2, f"  Dimensions ntax={len(taxa)};")
        lines.insert(3, "  Taxlabels " + " ".join(taxa) + ";")
        lines.insert(4, "End;")
    return "\n".join(lines)


def test_nexusreader_count_trees_simple():
    data = make_nexus(["(A,B);", "(C,D);"])
    f = io.StringIO(data)
    with NexusReader(f) as nf:
        assert nf.count_trees() == 2


def test_nexusreader_parse_simple():
    data = make_nexus(["(A,B);", "(B,A);"])
    f = io.StringIO(data)
    nf = NexusReader(f)
    trees = nf.parse()
    assert trees[0].newick() == "(A,B);"
    assert trees[1].newick() == "(B,A);"
    nf.close()


def test_nexusreader_next_and_has_next():
    data = make_nexus(["(A,B);", "(B,A);"])
    f = io.StringIO(data)
    nf = NexusReader(f)
    assert nf.has_next()
    t1 = nf.next()
    assert t1.newick() == "(A,B);"
    assert nf.has_next()
    t2 = nf.next()
    assert t2.newick() == "(B,A);"
    assert not nf.has_next()
    assert nf.next() is None
    nf.close()


def test_nexusreader_skip_next():
    data = make_nexus(["(A,B);", "(B,A);"])
    f = io.StringIO(data)
    nf = NexusReader(f)
    assert nf.has_next()
    nf.skip_next()
    assert nf.has_next()
    t = nf.next()
    assert t.newick() == "(B,A);"
    assert not nf.has_next()
    nf.close()


def test_nexusreader_empty():
    f = io.StringIO("")
    with NexusReader(f) as nf:
        assert nf.count_trees() == 0
        assert nf.parse() == []
        assert not nf.has_next()
        assert nf.next() is None


def test_nexusreader_translate_block():
    translate = (("1", "A"), ("2", "B"))
    trees = ["(1,2);"]
    data = make_nexus(trees, translate=translate)
    f = io.StringIO(data)
    with NexusReader(f) as nf:
        trees = nf.parse()
    assert len(trees) == 1
    assert set(nf.taxon_names) == {"A", "B"}
    assert trees[0].newick() == "(A,B);"


def test_nexusreader_with_taxlabels():
    taxa = ["A", "B", "C"]
    trees = ["(C,B,A);", "(B,C,A);"]
    data = make_nexus(trees, taxa=taxa)
    f = io.StringIO(data)
    nf = NexusReader(f, taxa)
    trees = nf.parse()
    nf.close()
    assert len(trees) == 2
    assert trees[0].taxon_names == taxa
    assert trees[1].taxon_names == taxa


def test_nexusreader_with_taxon_space():
    translate = [("1", "A A"), ("2", "B B"), ("3", "C C")]
    taxon_names = ["A A", "B B", "C C"]
    trees = ["(3,2,1);", "(2,3,1);"]
    data = make_nexus(trees, translate=translate)
    f = io.StringIO(data)
    nf = NexusReader(f, taxon_names)
    trees = nf.parse()
    nf.close()
    assert len(trees) == 2
    assert trees[0].taxon_names == taxon_names
    assert trees[1].taxon_names == taxon_names


def test_nexusreader_raise_on_different_taxa():
    data = "#NEXUS\nBegin trees;\ntree t1 = (A,B);\n  tree t2 = (C,A);\nEnd;"
    f = io.StringIO(data)
    nf = NexusReader(f)
    assert nf.count_trees() == 2
    with pytest.raises(Exception):
        nf.parse()
    nf.close()


def test_nexusreader_comments():
    data = "#NEXUS\nBegin trees;\n[coment]\ntree t1 = (A,B);\n  tree t2 = (B,A);\nEnd;"
    f = io.StringIO(data)
    nf = NexusReader(f)
    assert nf.count_trees() == 2
    trees = nf.parse()
    assert trees[0].newick() == "(A,B);"
    assert trees[1].newick() == "(B,A);"
    nf.close()


def test_nexusreader_skip_next_on_empty():
    f = io.StringIO("")
    nf = NexusReader(f)
    nf.skip_next()  # Should not raise
    nf.close()


def test_nexusreader_skip_next_multiple():
    data = make_nexus(["(A,B,C);", "(C,A,B);", "(B,A,C);"])
    f = io.StringIO(data)
    nf = NexusReader(f)
    nf.skip_next()
    nf.skip_next()
    assert nf.has_next()
    t = nf.next()
    assert t.newick() == "(B,A,C);"
    assert not nf.has_next()
    nf.close()


def test_nexusreader_with_comment():
    data = """#NEXUS
    Begin trees;
    tree t1 [&lnl=-1,prior=0.1,a=[0,2]] = [&R] (A,B);
    tree t2 [&lnl=-2,prior=0.2,a=[1,3]] = [&R] (B,A);
    End;
    """

    f = io.StringIO(data)
    nf = NexusReader(f)
    trees = nf.parse()
    nf.close()
    assert trees[0].newick() == "(A,B);"
    assert trees[1].newick() == "(B,A);"
    assert trees[0].comment == "[&lnl=-1,prior=0.1,a=[0,2]]"
    assert trees[1].comment == "[&lnl=-2,prior=0.2,a=[1,3]]"


def test_nexuswriter_save_and_read(tmp_path):
    # Create a simple tree and save it
    tree = Tree.from_newick("(A,B);", ["A", "B"])
    file_path = tmp_path / "test.nex"
    NexusWriter.save(str(file_path), tree, include_translate=False)
    # Read back and check content
    with open(file_path) as f:
        content = f.read()
    assert "#NEXUS" in content
    assert "Begin trees;" in content
    assert "Translate" not in content
    assert "tree STATE_1" in content
    assert "(A,B);" in content
    assert "End;" in content


def test_nexuswriter_save_multiple_trees(tmp_path):
    trees = [
        Tree.from_newick("(A,B);", ["A", "B"]),
        Tree.from_newick("(B,A);", ["A", "B"]),
    ]
    file_path = tmp_path / "multi.nex"
    NexusWriter.save(str(file_path), trees, tree_prefix="T_")
    with open(file_path) as f:
        content = f.read()
    expected = textwrap.dedent(
        """\
        #NEXUS
        Begin trees;
        Translate
          1 A,
          2 B;
        tree T_1 = [&R] (1,2);
        tree T_2 = [&R] (2,1);
        End;
    """
    )
    assert expected == content


def test_nexuswriter_include_tree_comment(tmp_path):
    tree = Tree.from_newick("(A,B);")
    tree.comment = "[&mycomment]"
    file_path = tmp_path / "comment.nex"
    NexusWriter.save(
        str(file_path), tree, include_tree_comment=True, include_translate=False
    )
    with open(file_path) as f:
        content = f.read()

    assert "tree STATE_1 [&mycomment] = [&R] (A,B);" in content


def test_nexuswriter_context_manager(tmp_path):
    tree = Tree.from_newick("(A,B);")
    file_path = tmp_path / "ctx.nex"
    with NexusWriter(str(file_path), include_translate=False) as writer:
        writer.begin_block("trees")
        writer.write([tree])
        writer.end_block()
    with open(file_path) as f:
        content = f.read()
    assert "tree STATE_1 = [&R] (A,B);" in content


def test_nexuswriter_tree_prefix_and_options(tmp_path):
    tree = Tree.from_newick("(A:0.1,B:0.2);")
    file_path = tmp_path / "prefix.nex"
    options = {
        "tree_prefix": "MYTREE_",
        "include_translate": True,
        "include_branch_lengths": False,
    }
    NexusWriter.save(
        str(file_path),
        tree,
        **options,
    )
    with open(file_path) as f:
        content = f.read()
    assert "tree MYTREE_1 = [&R] (1,2);" in content
