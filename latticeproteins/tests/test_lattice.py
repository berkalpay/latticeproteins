from unittest import TestCase
import pickle

from latticeproteins import Lattice, Protein


def test_num_conformations():
    assert len(list(Lattice(L=6).ensemble.conformations)) == 36


class TestFolding(TestCase):
    def setUp(self):
        self.lattice = Lattice(L=18)
        self.protein = Protein("CDEFFKKHCIERMFMCYW")
        self.lattice.fold(self.protein)

    def test_does_fold(self):
        self.assertIsNotNone(self.protein.native_state)

    def test_conformation(self):
        self.assertEqual(self.protein.native_state.bond_dirs, "URDRURDDDDLUULLDR")

    def test_stability(self):
        self.assertEqual(round(self.protein.stability(temp=1.0), 2), -1.04)


def test_lattices_equal():
    assert Lattice(L=10) == Lattice(L=10)


def test_lattices_unequal():
    assert Lattice(L=9) != Lattice(L=10)
