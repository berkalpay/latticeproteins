from unittest import TestCase
import pickle

from latticeproteins import Lattice, Protein


def test_num_conformations():
    assert len(list(Lattice(L=6).ensemble.conformations())) == 36


class TestFolding(TestCase):
    def setUp(self):
        #self.lattice = Lattice(L=18)
        self.lattice = pickle.load(open("../lattice18.pickle", "rb"))  # TODO: standardize path
        self.protein = Protein("CDEFFKKHCIERMFMCYW")
        self.lattice.fold(self.protein)

    def test_does_fold(self):
        assert self.protein.native_state

    def test_conformation(self):
        assert self.protein.native_state.bond_dirs == "URDRURDDDDLUULLDR"

    def test_stability(self):
        assert round(self.protein.stability(temp=1.0), 2) == -1.04
