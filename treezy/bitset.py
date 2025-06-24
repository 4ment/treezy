# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT


from typing import List, Union


class BitSet:
    """A simple bitset supporting bitwise and bit manipulation operations.

    Attributes
    ----------
    size
        The number of bits in the bitset.
    value
        The integer value representing the bitset, where each bit corresponds
        to a position in the set.

    Example
    -------
    >>> a = BitSet(4)
    >>> a.to_list()
    [0, 0, 0, 0]
    >>> a[1] = True
    >>> a.to_list()
    [0, 1, 0, 0]
    >>> print(a)
    0100
    >>> b = BitSet(4, 0b1010)
    >>> print(b)
    0101
    >>> a |= b
    >>> print(a)
    0101
    """

    size: int
    value: int

    def __init__(self, size: int = 0, value: int = 0):
        """Initialize a BitSet with a given size and value.

        Parameters
        ----------
        size
            The number of bits in the bitset.
        value
            The initial value of the bitset, defaults to 0.
        """
        if size <= 0:
            raise ValueError("Size must be a positive integer.")
        self.size = size
        self.value = value & ((1 << size) - 1) if size > 0 else value

    def __repr__(self):
        return f"BitSet(size={self.size}, value={bin(self.value)})"

    def __str__(self) -> str:
        return ''.join(str((self.value >> i) & 1) for i in range(self.size))

    def __and__(self, other):
        size = max(self.size, other.size)
        return BitSet(size, self.value & other.value)

    def __or__(self, other):
        size = max(self.size, other.size)
        return BitSet(size, self.value | other.value)

    def __xor__(self, other):
        size = max(self.size, other.size)
        return BitSet(size, self.value ^ other.value)

    def __invert__(self):
        mask = (1 << self.size) - 1
        return BitSet(self.size, ~self.value & mask)

    def __lshift__(self, n):
        return BitSet(self.size, (self.value << n) & ((1 << self.size) - 1))

    def __rshift__(self, n):
        return BitSet(self.size, self.value >> n)

    def __eq__(self, other):
        return self.size == other.size and self.value == other.value

    def __hash__(self):
        return hash((self.size, self.value))

    def __len__(self) -> int:
        return self.size

    def __contains__(self, pos: int) -> bool:
        return self.test(pos)

    def set(self, pos: int) -> None:
        """Set the bit at the given position.

        Parameters
        ----------
        pos
            The position of the bit to set.

        Raises
        ------
        IndexError
            If the position is out of range.
        """
        if pos < self.size:
            self.value |= 1 << pos
        else:
            raise IndexError(f"Position out of range ({pos} >= {self.size})")

    def clear(self, pos: int) -> None:
        """Clear the bit at the given position.

        Parameters
        ----------
        pos
            the position of the bit to clear.

        Raises
        ------
        IndexError
            If the position is out of range.
        """
        if pos < self.size:
            self.value &= ~(1 << pos)
        else:
            raise IndexError(f"Position out of range ({pos} >= {self.size})")

    def flip(self, pos: int = None) -> None:
        """Flip the bit at the given position or all bits if pos is None.

        Parameters
        ----------
        pos
            The position of the bit to flip. If None, all bits are flipped.

        Raises
        ------
        IndexError
            If the position is out of range.
        """
        if pos is None:
            self.value = (1 << self.size) - 1 - self.value
        elif pos < self.size:
            self.value ^= 1 << pos
        else:
            raise IndexError(f"Position out of range ({pos} >= {self.size})")

    def test(self, pos: int) -> bool:
        """Test if the bit at the given position is set.

        Parameters
        ----------
        pos
            The position of the bit to test.

        Returns
        -------
        bool
            True if the bit is set, False otherwise.

        Raises
        ------
            IndexError: If the position is out of range.
        """
        if pos < self.size:
            return (self.value & (1 << pos)) != 0
        else:
            raise IndexError(f"Position out of range ({pos} >= {self.size})")

    def count(self) -> int:
        """Count the number of bits set to 1.

        Returns
        -------
        int
            The number of bits set to 1.
        """
        return bin(self.value).count('1')

    def to_list(self) -> List[int]:
        """Convert the bitset to a list of bits.

        Returns
        -------
        List[int]
            A list of bits (0 or 1) representing the bitset.
        """
        return [(self.value >> i) & 1 for i in range(self.size)]

    def __getitem__(self, pos: int) -> bool:
        """Get the value of the bit at the given position.

        Parameters
        ----------
        pos
            The position of the bit to get.

        Returns
        -------
        bool
            True if the bit is set, False otherwise.

        Raises
        ------
        IndexError
            If the position is out of range.
        """
        return self.test(pos)

    def __setitem__(self, pos: int, val: Union[int, bool]) -> None:
        """Set the value of the bit at the given position.

        Parameters
        ----------
        pos
            The position of the bit to set.
        val
            The value to set (True or False).

        Raises
        ------
        IndexError
            If the position is out of range.
        """
        if val:
            self.set(pos)
        else:
            self.clear(pos)
