# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

from abc import ABC, abstractmethod
from typing import Union

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

    def __init__(self, weighted: bool = False):
        """Initializes a RobinsonFouldsMetric.

        Parameters
        ----------
        weighted
            If True, the distance is weighted by the branch lengths of the splits.
        """
        super().__init__()
        self._weighted = weighted

    def compute_from_splits(
        self, splits1: Union[dict, set], splits2: Union[dict, set]
    ) -> float:
        """Compute the Robinson-Foulds distance from two collections of splits.

        If the splits are sets, they are treated as unweighted.
        If they are dictionaries, the values are treated as weights for the splits.

        Parameters
        ----------
        splits1
            splits from the first tree.
        splits2
            splits from the second tree.

        Returns
        -------
        float
            The Robinson-Foulds distance between the two collections of splits.
        """
        if isinstance(splits1, set):
            return float(len(splits1 ^ splits2))
        elif isinstance(splits1, dict):
            return float(
                sum(
                    abs(splits1.get(s, 0.0) - splits2.get(s, 0.0))
                    for s in set(splits1.keys()).union(splits2.keys())
                )
            )
        else:
            raise ValueError("splits1 and splits2 must be either sets or dictionaries.")

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
        if tree1.is_rooted != tree2.is_rooted:
            raise ValueError("Both trees must be either rooted or unrooted.")

        if tree1.is_rooted:
            if not self._weighted:
                splits1 = {
                    n.descendant_bitset.value
                    for n in tree1.nodes
                    if not n.is_leaf and not n.is_root
                }
                splits2 = {
                    n.descendant_bitset.value
                    for n in tree2.nodes
                    if not n.is_leaf and not n.is_root
                }
            else:
                splits1 = {
                    n.descendant_bitset.value: n.distance
                    for n in tree1.nodes
                    if not n.is_leaf and not n.is_root
                }
                splits2 = {
                    n.descendant_bitset.value: n.distance
                    for n in tree2.nodes
                    if not n.is_leaf and not n.is_root
                }
        else:
            if not self._weighted:
                splits1 = {
                    frozenset((n.descendant_bitset.value, (~n.descendant_bitset).value))
                    for n in tree1.nodes
                    if not n.is_leaf and not n.is_root
                }
                splits2 = {
                    frozenset((n.descendant_bitset.value, (~n.descendant_bitset).value))
                    for n in tree2.nodes
                    if not n.is_leaf and not n.is_root
                }
            else:
                splits1 = {
                    frozenset(
                        (n.descendant_bitset.value, (~n.descendant_bitset).value)
                    ): n.distance
                    for n in tree1.nodes
                    if not n.is_leaf and not n.is_root
                }
                splits2 = {
                    frozenset(
                        (n.descendant_bitset.value, (~n.descendant_bitset).value)
                    ): n.distance
                    for n in tree2.nodes
                    if not n.is_leaf and not n.is_root
                }
        return self.compute_from_splits(splits1, splits2)
