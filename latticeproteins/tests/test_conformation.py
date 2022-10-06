from unittest import TestCase

from latticeproteins import Conformation, Lattice


class ConformationTestCase(TestCase):
    def setUp(self):
        self.conformation = Conformation("UUURURDDDL")
        self.true_locs = {(0, 0), (0, 1), (0, 2), (0, 3), (1, 3), (1, 4),
                          (2, 4), (2, 3), (2, 2), (2, 1), (1, 1)}

    def test_width(self):
        self.assertEqual(self.conformation.width, 3)

    def test_height(self):
        self.assertEqual(self.conformation.height, 5)

    def test_not_overlapping(self):
        self.assertFalse(self.conformation.overlapping)

    def test_locations(self):
        self.assertEqual(set(self.conformation.locations), self.true_locs)

    def test_shifted_locations(self):
        shifted_conformation = Conformation(self.conformation.bond_dirs)
        shifted_conformation.location_delta = (-2, 1)
        true_delta_locs = set([(x - 2, y + 1) for x, y in self.true_locs])
        self.assertEqual(set(shifted_conformation.locations), true_delta_locs)

    def test_rotating_zero_times(self):
        self.assertEqual(self.conformation.rotated_clockwise(0).bond_dirs, self.conformation.bond_dirs)

    def test_rotating(self):
        self.assertEqual(self.conformation.rotated_clockwise(3).bond_dirs, "LLLULURRRD")


def test_rotating():
    assert Conformation("UUU").rotated_clockwise(3).bond_dirs == "LLL"


def test_overlapping():
    assert Conformation("UUURURDDDLL").overlapping


def test_non_overlapping():
    assert not Conformation("UUURURUUR").overlapping


def test_num_contacts():
    assert len(Conformation("URDRUULLLDDDRRDLL").contacts) == 10


def test_max_num_contacts():
    assert max(len(conf.contacts) for conf in Lattice(L=6).ensemble.conformations) == 2
