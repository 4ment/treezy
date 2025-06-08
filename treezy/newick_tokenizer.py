# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

from typing import Generator, NamedTuple


class Token(NamedTuple):
    """A token representing a component of a Newick string."""

    type: str
    value: str


def tokenize_newick(newick: str) -> Generator[Token, None, None]:
    """Tokenize a Newick string into a sequence of tokens.

    This function processes a Newick-formatted string and yields tokens
    representing different components of the Newick format, such as
    parentheses, commas, colons, semicolons, comments, names, and branch lengths.

    Rules for a valid Newick string:
    - taxon and internal node names cannot contain special characters like ``(),:;[]``.
    - comments are enclosed in square brackets ``[]`` and can contain any characters
      except for the closing bracket ``]``.
    - branch lengths can be integers or floating-point numbers, including scientific
      notation.
    - a semicolon ``;`` is used to terminate the Newick string.
    - there should be an equal number of ``(`` and ``)``.
    - comments placed before a taxon name or closing bracket represent node comments.
    - comments placed between a colon ``:`` and the branch length represent
      branch comments.

    The function yields tokens in the order they appear in the Newick string.

    The tokens yielded are:
        - ``LPAREN``: for '('
        - ``RPAREN``: for ')'
        - ``COMMA``: for ','
        - ``COLON``: for ':'
        - ``SEMICOLON``: for ';'
        - ``COMMENT``: for comments enclosed in square brackets '[]'
        - ``NAME``: for taxon names or internal node names
        - ``NUMBER``: for branch lengths

    This function raises a ValueError if it encounters an unclosed comment.

    Parameters
    ----------
    newick
        The Newick string to tokenize.

    Raises
    ------
    ValueError
        If an unclosed comment is found in the Newick string.

    Returns
    -------
    Generator[Token]
        A generator yielding tokens parsed from the Newick string.
    """
    i = 0
    n = len(newick)
    expect_a_number = False
    opened_parentheses = 0

    while i < n:
        c = newick[i]

        if c == '(':
            opened_parentheses += 1
            yield Token('LPAREN', c)
        elif c == ')':
            opened_parentheses -= 1
            if opened_parentheses < 0:
                raise ValueError("Unmatched closing parenthesis in Newick string")
            yield Token('RPAREN', c)
        elif c == ',':
            yield Token('COMMA', c)
        elif c == ':':
            yield Token('COLON', c)
            expect_a_number = True
        elif c == ';':
            yield Token('SEMICOLON', c)
        elif c == '[':
            start = i
            i += 1
            count = 1
            while i < n and count > 0:
                if newick[i] == ']':
                    count -= 1
                elif newick[i] == '[':
                    count += 1
                i += 1
            if i >= n:
                raise ValueError(f"Unclosed comment at position {start}")
            yield Token('COMMENT', newick[start:i])
            i -= 1
        else:
            start = i
            while i < n and newick[i] not in '():;,[':
                i += 1
            value = newick[start:i]
            i -= 1

            if expect_a_number:
                yield Token('NUMBER', value)
                expect_a_number = False
            else:
                yield Token('NAME', value)

        i += 1
