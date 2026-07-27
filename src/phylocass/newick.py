"""A small Newick reader for rooted phylogenetic trees.

Supports leaf and internal labels, quoted labels, branch lengths, bootstrap
values and ``[...]`` comments.  Branch lengths are parsed and discarded --
Cass only cares about topology.
"""

from __future__ import annotations

from typing import Iterable

from .network import Network

__all__ = ["parse_newick", "parse_newick_file", "read_trees"]


class NewickError(ValueError):
    """Raised when a Newick string cannot be parsed."""


class _Parser:
    def __init__(self, text: str) -> None:
        self.text = text
        self.pos = 0

    def error(self, msg: str) -> NewickError:
        return NewickError(f"{msg} at position {self.pos} in {self.text[:60]!r}")

    def skip(self) -> None:
        while self.pos < len(self.text):
            ch = self.text[self.pos]
            if ch.isspace():
                self.pos += 1
            elif ch == "[":
                depth = 0
                while self.pos < len(self.text):
                    if self.text[self.pos] == "[":
                        depth += 1
                    elif self.text[self.pos] == "]":
                        depth -= 1
                        if depth == 0:
                            self.pos += 1
                            break
                    self.pos += 1
                else:
                    raise self.error("unterminated comment")
            else:
                return

    def peek(self) -> str:
        self.skip()
        return self.text[self.pos] if self.pos < len(self.text) else ""

    def read_label(self) -> str:
        self.skip()
        if self.pos < len(self.text) and self.text[self.pos] in "'\"":
            quote = self.text[self.pos]
            self.pos += 1
            out = []
            while self.pos < len(self.text):
                ch = self.text[self.pos]
                if ch == quote:
                    self.pos += 1
                    if self.pos < len(self.text) and self.text[self.pos] == quote:
                        out.append(quote)
                        self.pos += 1
                        continue
                    return "".join(out)
                out.append(ch)
                self.pos += 1
            raise self.error("unterminated quoted label")
        start = self.pos
        while self.pos < len(self.text) and self.text[self.pos] not in "(),:;[]":
            self.pos += 1
        return self.text[start:self.pos].strip().replace("_", " ")

    def skip_length(self) -> None:
        if self.peek() == ":":
            self.pos += 1
            self.skip()
            while self.pos < len(self.text) and self.text[self.pos] not in "(),;[]":
                self.pos += 1

    def parse(self, tree: Network) -> int:
        self.skip()
        if self.peek() == "(":
            self.pos += 1
            node = tree.new_node()
            while True:
                child = self.parse(tree)
                tree.add_edge(node, child)
                nxt = self.peek()
                if nxt == ",":
                    self.pos += 1
                    continue
                if nxt == ")":
                    self.pos += 1
                    break
                raise self.error("expected ',' or ')'")
            self.read_label()  # internal labels / support values are ignored
            self.skip_length()
            return node
        name = self.read_label()
        if not name:
            raise self.error("expected a leaf label")
        self.skip_length()
        return tree.new_node(label=frozenset({name}))


def parse_newick(text: str) -> Network:
    """Parse a single Newick string into a :class:`~phylocass.network.Network`."""
    parser = _Parser(text)
    tree = Network()
    root = parser.parse(tree)
    if parser.peek() == ";":
        parser.pos += 1
    parser.skip()
    if parser.pos < len(parser.text):
        raise parser.error("trailing characters after tree")
    if len(tree.roots()) != 1 or tree.roots()[0] != root:
        raise NewickError("parsed tree does not have a unique root")
    seen: set[frozenset] = set()
    for v in tree.leaves():
        lbl = tree.label.get(v)
        if lbl is None:
            raise NewickError("tree has an unlabelled leaf")
        if lbl in seen:
            raise NewickError(f"duplicate leaf label {sorted(lbl)[0]!r}")
        seen.add(lbl)
    return tree


def read_trees(text: str) -> list[Network]:
    """Parse every ``;``-terminated Newick tree in ``text``.

    Blank lines and lines starting with ``#`` are ignored.
    """
    trees: list[Network] = []
    buf: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        buf.append(line)
        joined = " ".join(buf)
        while ";" in joined:
            head, joined = joined.split(";", 1)
            trees.append(parse_newick(head + ";"))
        buf = [joined] if joined.strip() else []
    if buf and " ".join(buf).strip():
        trees.append(parse_newick(" ".join(buf)))
    return trees


def parse_newick_file(path: str) -> list[Network]:
    with open(path, "r", encoding="utf-8") as fh:
        return read_trees(fh.read())
