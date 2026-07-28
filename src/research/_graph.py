"""Small graph utilities backing the graphic and transversal constructors.

A single union-find serves both the forest test behind graphic matroids and
the circuit-connectivity components of positroids; the bipartite matching is
the independence oracle for transversal matroids.
"""

from collections.abc import Hashable, Iterator, Sequence

from research._bitmask import bits

__all__ = [
    "UnionFind",
    "has_matching",
    "is_forest",
]


class UnionFind[V: Hashable]:
    """A disjoint-set forest over hashable items, with path compression.

    Items register on first touch; :meth:`groups` reports the partition of
    everything registered so far.
    """

    def __init__(self) -> None:
        """Start with no items registered."""
        self._parent: dict[V, V] = {}

    def add(self, item: V) -> None:
        """Register ``item`` as a singleton component if it is unseen."""
        self._parent.setdefault(item, item)

    def find(self, item: V) -> V:
        """Return the representative of ``item``'s component, registering it."""
        self.add(item)
        root = item
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[item] != root:
            self._parent[item], item = root, self._parent[item]
        return root

    def union(self, a: V, b: V) -> bool:
        """Merge the components of ``a`` and ``b``.

        Returns:
            Whether the two were in distinct components before the merge.
        """
        root_a, root_b = self.find(a), self.find(b)
        if root_a == root_b:
            return False
        self._parent[root_a] = root_b
        return True

    def groups(self) -> Iterator[frozenset[V]]:
        """Yield the current components over every registered item.

        The iterator is single-use and reflects the partition at call time.
        """
        by_root: dict[V, set[V]] = {}
        for item in list(self._parent):
            by_root.setdefault(self.find(item), set()).add(item)
        for group in by_root.values():
            yield frozenset(group)


def is_forest[V: Hashable](edges: Sequence[tuple[V, V]]) -> bool:
    """Return whether the given edge multiset is acyclic (a forest)."""
    dsu = UnionFind[V]()
    return all(dsu.union(u, v) for u, v in edges)


def has_matching(mask: int, set_masks: Sequence[int]) -> bool:
    """Return whether every element bit of ``mask`` gets a distinct set.

    Kuhn's augmenting-path bipartite matching between the elements of
    ``mask`` and the members of ``set_masks`` containing them.
    """
    matched: dict[int, int] = {}

    def augment(element: int, visited: set[int]) -> bool:
        for set_index, set_mask in enumerate(set_masks):
            if set_mask >> element & 1 and set_index not in visited:
                visited.add(set_index)
                if set_index not in matched or augment(matched[set_index], visited):
                    matched[set_index] = element
                    return True
        return False

    return all(augment(element, set()) for element in bits(mask))
