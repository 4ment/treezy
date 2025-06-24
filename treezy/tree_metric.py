# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

from abc import ABC, abstractmethod
from typing import Set

from treezy.tree import Tree


class TreeMetric(ABC):
    """Abstract base class for tree metrics."""

    @abstractmethod
    def compute(self, tree1: Tree, tree2: Tree) -> float:
        """Compute the metric between two trees."""


class RobinsonFouldsMetric(TreeMetric):
    r"""Robinson-Foulds metric for comparing two trees.

    This metric counts the number of splits that are present in one tree but not
    in the other.

    Let :math:`\mathcal{B}(T_1)` be the set of splits in tree :math:`T_1` and
    :math:`\mathcal{B}(T_2)` be the set of splits in tree :math:`T_2`.
    The Robinson-Foulds distance is defined as:

    .. math::

       \mathrm{RF}(T_1, T_2) = |\mathcal{B}(T_1) \setminus \mathcal{B}(T_2)| +
       |\mathcal{B}(T_2) \setminus \mathcal{B}(T_1)|

    or equivalently:

    .. math::

       \mathrm{RF}(T_1, T_2) = 2 \cdot | \mathcal{B}(T_1) \cup \mathcal{B}(T_2)| -
       2 \cdot | \mathcal{B}(T_1) \cap \mathcal{B}(T_2) |.
    """

    def compute_from_sets(self, set1: Set, set2: Set) -> float:
        """Compute the Robinson-Foulds distance from two sets of splits.


        Parameters
        ----------
        set1
            Set of splits from the first tree.
        set2
            Set of splits from the second tree.

        Returns
        -------
        float
            The Robinson-Foulds distance between the two sets.
        """
        shared = len(set1 & set2)
        total = len(set1) + len(set2)
        return float(total - 2 * shared)

    def compute(self, tree1: Tree, tree2: Tree) -> float:
        """Compute the Robinson-Foulds distance between two trees.

        It assumes that the descentant bitsets of the nodes of the trees
        are already computed and available in the nodes:

            .. code-block:: python

                tree1.compute_descendant_bitset()
                tree2.compute_descendant_bitset()

        Parameters
        ----------
        tree1
            The first tree.
        tree2
            The second tree.
        Returns
        -------
        float
            The Robinson-Foulds distance between the two trees.
        """
        set1 = {
            n.descendant_bitset for n in tree1.nodes if not n.is_leaf and not n.is_root
        }
        set2 = {
            n.descendant_bitset for n in tree2.nodes if not n.is_leaf and not n.is_root
        }
        return self.compute_from_sets(set1, set2)
