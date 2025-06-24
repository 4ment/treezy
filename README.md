# treezy

[![Tests](https://github.com/4ment/treezy/actions/workflows/test.yml/badge.svg)](https://github.com/4ment/treezy/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
![PyPI](https://img.shields.io/pypi/v/treezy)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/treezy)

treezy is a simple Python library for manipulating phylogenetic trees, including reading from and writing to files.

## Features

- Parse and manipulate phylogenetic trees
- Read and write trees with branch annotations and comments
- Lightweight and easy to integrate into existing Python projects

## Installation

```bash
pip install treezy
```

To build treezy from source you can run
```bash
git clone https://github.com/4ment/treezy
pip install treezy/
```

## Usage

```python
from treezy.tree import Tree

newick = '((A:[&rate=0.1,cat=0]1,B:[&rate=0.1,cat=0]2):[&rate=0.2,cat=1]3,C:[&rate=0.2,cat=1]4);'
tree = Tree.from_newick(newick)
tree.parse_branch_comment({'rate': lambda rate: float(rate)})
for node in tree:
    if not node.is_root:
        node.distance *= node.branch_annotations['rate']
print(tree)
# ((A:0.1,B:0.2):0.6,C:0.8);
```

## Contributing

Contributions are welcome! Please open issues or submit pull requests.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
