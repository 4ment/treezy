# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

class BitSet:
    def __init__(self, size=0, value=0):
        self.size = size
        self.value = value & ((1 << size) - 1) if size > 0 else value

    def __repr__(self):
        return f"BitSet(size={self.size}, value={bin(self.value)})"

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

    def set(self, pos):
        if pos < self.size:
            self.value |= 1 << pos

    def clear(self, pos):
        if pos < self.size:
            self.value &= ~(1 << pos)

    def toggle(self, pos: int) -> None:
        if pos < self.size:
            self.value ^= 1 << pos

    def test(self, pos):
        if pos < self.size:
            return (self.value & (1 << pos)) != 0
        return False

    def count(self):
        return bin(self.value).count('1')

    def to_list(self):
        return [(self.value >> i) & 1 for i in range(self.size)]

    def __getitem__(self, pos):
        return self.test(pos)

    def __setitem__(self, pos, val):
        if val:
            self.set(pos)
        else:
            self.clear(pos)
