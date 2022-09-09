from unittest import TestCase
from latticeproteins.core import Lattice, Protein, Conformation
from latticeproteins.thermodynamics import stability


def test_num_conformations():
    assert len(Lattice(L=6).conformations) == 36


class TestFolding(TestCase):
    def setUp(self):
        self.lattice = Lattice(L=18)
        self.protein = Protein("CDEFFKKHCIERMFMCYW")
        self.lattice.fold(self.protein)

    def test_does_fold(self):
        conf = self.protein.conformation
        assert conf and isinstance(conf, Conformation)

    def test_conformation(self):
        assert self.protein.conformation == Conformation("URDDLLDRRRRULURDL")

    def test_stability(self):
        assert stability(self.lattice.energy(self.protein.seq, self.protein.conformation),
                         sum(self.lattice.all_energies(self.protein.seq)),
                         1.0) == -1.04
