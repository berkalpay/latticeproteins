from unittest import TestCase

from latticeproteins import Lattice, Protein, Conformation


def test_binding_energy():
    protein = Protein("KIR", conformations=[Conformation("UU")], lattice=Lattice(L=3))
    ligand = Protein("NI", conformations=[Conformation("U")], lattice=Lattice(L=2))
    ligand.native_state = ligand.native_state.rotated_clockwise(2)
    ligand.native_state.location_delta = (-1, 2)
    print(ligand.native_state.locations)
    contact_set = protein.native_state.contacts_with(ligand.native_state)
    total_contact_energy = Lattice(L=3).sum_contact_energy(protein.seq, contact_set, ligand.seq)  # TODO: shouldn't have to do this
    assert round(total_contact_energy, 2) == -7.63


class TestBindingSmall(TestCase):
    def setUp(self):
        protein = Protein("KIR", conformations=[Conformation("UU")], lattice=Lattice(L=3))
        ligand = Protein("NI", conformations=[Conformation("U")], lattice=Lattice(L=2))
        self.binding_energy, self.positioning_info = protein.bind(ligand)

    def test_min_binding_energy(self):
        self.assertEqual(round(self.binding_energy, 2), -7.63)

    def test_min_binding_positioning(self):
        for rotations, location_delta in self.positioning_info:
            self.assertEqual(rotations, 2)
            self.assertIn(location_delta, [(-1, 2), (1, 2)])


class TestBindingLarge(TestCase):
    def setUp(self):
        protein = Protein("CDEFFKKHCIERMFMCYW", conformations=[Conformation("URDRURDDDDLUULLDR")], lattice=Lattice(L=18))
        ligand = Protein("HDGEKGKA", conformations=[Conformation("URUUULL")], lattice=Lattice(L=2))
        self.binding_energy, self.positioning_info = protein.bind(ligand)

    def test_min_binding_energy(self):
        self.assertEqual(round(self.binding_energy, 2), -17.90)

    def test_min_binding_positioning(self):
        for rotations, location_delta in self.positioning_info:
            self.assertEqual(rotations, 1)
            self.assertEqual(location_delta, (0, -3))
