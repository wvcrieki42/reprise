"""Disease-ontology helpers for novelty radius computation."""
from __future__ import annotations
from collections import defaultdict
import pandas as pd


class DiseaseOntology:
    """Parent/child navigation over an EFO-style ontology edge list."""

    def __init__(self, ontology: pd.DataFrame):
        self.parents: dict[str, set[str]] = defaultdict(set)
        self.children: dict[str, set[str]] = defaultdict(set)
        self.name: dict[str, str] = {}
        for efo, dname, parent in ontology[["efo_id", "disease_name", "parent_efo_id"]].itertuples(index=False):
            self.name.setdefault(efo, dname)
            if parent:
                self.parents[efo].add(parent)
                self.children[parent].add(efo)

    def neighborhood(self, efo_ids: set[str], radius: int) -> dict[str, int]:
        """Return {efo_id: hop_distance} for all nodes within `radius` of the seeds.

        Distance 0 = the seed itself; 1 = direct parent or child; etc.
        Used to mark a disease as 'known' if it is close to an existing indication.
        """
        dist: dict[str, int] = {e: 0 for e in efo_ids}
        frontier = set(efo_ids)
        for hop in range(1, radius + 1):
            nxt: set[str] = set()
            for node in frontier:
                for nb in self.parents.get(node, set()) | self.children.get(node, set()):
                    if nb not in dist:
                        dist[nb] = hop
                        nxt.add(nb)
            frontier = nxt
            if not frontier:
                break
        return dist
