# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

import pytest

from treezy.bitset import BitSet


def test_init_and_repr():
    b = BitSet(4, 5)
    assert b.size == 4
    assert b.value == 5
    assert repr(b) == "BitSet(size=4, value=0b101)"


def test_and_or_xor():
    a = BitSet(4, 0b1010)
    b = BitSet(4, 0b1100)
    assert (a & b) == BitSet(4, 0b1000)
    assert (a | b) == BitSet(4, 0b1110)
    assert (a ^ b) == BitSet(4, 0b0110)


def test_invert():
    a = BitSet(4, 0b1010)
    assert ~a == BitSet(4, 0b0101)


def test_lshift_rshift():
    a = BitSet(4, 0b0011)
    assert (a << 1) == BitSet(4, 0b0110)
    assert (a << 2) == BitSet(4, 0b1100)
    assert (a << 3) == BitSet(4, 0b1000)
    assert (a << 4) == BitSet(4, 0b0000)
    assert (a >> 1) == BitSet(4, 0b0001)
    assert (a >> 2) == BitSet(4, 0b0000)


def test_eq():
    a = BitSet(4, 0b1010)
    b = BitSet(4, 0b1010)
    c = BitSet(4, 0b0101)
    d = BitSet(3, 0b010)
    assert a == b
    assert a != c
    assert a != d


def test_set_clear_flip_test():
    b = BitSet(4, 0)
    b.set(1)
    assert b.value == 0b0010
    b.set(3)
    assert b.value == 0b1010
    b.clear(1)
    assert b.value == 0b1000
    b.flip(0)
    assert b.value == 0b1001
    b.flip(3)
    assert b.value == 0b0001
    assert b.test(0) is True
    assert b.test(1) is False
    with pytest.raises(IndexError):
        assert b.test(4) is False  # out of range


def test_count_to_list():
    b = BitSet(5, 0b10101)
    assert b.count() == 3
    assert b.to_list() == [1, 0, 1, 0, 1]


def test_getitem_setitem():
    b = BitSet(4, 0)
    b[2] = 1
    assert b[2] is True
    b[2] = 0
    assert b[2] is False
    b[3] = True
    assert b[3] is True
    b[3] = False
    assert b[3] is False


def test_size_zero():
    with pytest.raises(ValueError):
        BitSet(0)
