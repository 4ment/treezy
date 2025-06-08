# Copyright 2025 Mathieu Fourment
# SPDX-License-Identifier: MIT

from treezy.newick_tokenizer import Token, tokenize_newick


def test_tokenize_newick_with_comments():
    newick = "((A[&key1={1,2}]:-0.1,1:[&key2=[1,2]]2)[&key3=3]:1e-2[&key4=val],C:1E+2);"
    tokens = tokenize_newick(newick)
    expected = [
        Token(type='LPAREN', value='('),
        Token(type='LPAREN', value='('),
        Token(type='NAME', value='A'),
        Token(type='COMMENT', value='[&key1={1,2}]'),
        Token(type='COLON', value=':'),
        Token(type='NUMBER', value='-0.1'),
        Token(type='COMMA', value=','),
        Token(type='NAME', value='1'),
        Token(type='COLON', value=':'),
        Token(type='COMMENT', value='[&key2=[1,2]]'),
        Token(type='NUMBER', value='2'),
        Token(type='RPAREN', value=')'),
        Token(type='COMMENT', value='[&key3=3]'),
        Token(type='COLON', value=':'),
        Token(type='NUMBER', value='1e-2'),
        Token(type='COMMENT', value='[&key4=val]'),
        Token(type='COMMA', value=','),
        Token(type='NAME', value='C'),
        Token(type='COLON', value=':'),
        Token(type='NUMBER', value='1E+2'),
        Token(type='RPAREN', value=')'),
        Token(type='SEMICOLON', value=';'),
    ]
    for idx, token in enumerate(tokens):
        assert (
            token == expected[idx]
        ), f"Token mismatch at index {idx}: expected {expected[idx]}, got {token}"
