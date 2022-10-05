from functools import cached_property, lru_cache
from dataclasses import dataclass
from typing import List
from itertools import chain
from math import inf
from copy import deepcopy

import numpy as np
from scipy.special import logsumexp

from latticeproteins.interactions import miyazawa_jernigan


def next_monomer_location(location, bond_dir):
    bond_dir_to_dx = {'U': 0, 'R': 1, 'D': 0, 'L': -1}
    bond_dir_to_dy = {'U': 1, 'R': 0, 'D': -1, 'L': 0}
    return location[0] + bond_dir_to_dx[bond_dir], location[1] + bond_dir_to_dy[bond_dir]


def bond_dir_rotated_clockwise(bond_dir):
    return {'U': 'R', 'R': 'D', 'D': 'L', 'L': 'U'}[bond_dir]


class Conformation:
    def __init__(self, bond_dirs, location_delta=(0, 0)):
        self._bond_dirs = bond_dirs
        self._location_delta = location_delta

    @property
    def bond_dirs(self):
        return self._bond_dirs

    @bond_dirs.setter
    def bond_dirs(self, value):
        self._bond_dirs = value
        self.__dict__.pop("contacts", None)

    @property
    def location_delta(self):
        return self._location_delta

    @location_delta.setter
    def location_delta(self, value):
        self._location_delta = value

    def __getitem__(self, item):
        return self.bond_dirs[item]

    def __setitem__(self, key, value):
        self.bond_dirs = self.bond_dirs[:key] + value + self.bond_dirs[key+1:]
        self.__dict__.pop("contacts", None)

    def __iter__(self):
        yield from self.bond_dirs

    def __eq__(self, other):
        return self.bond_dirs == other.bond_dirs

    def __hash__(self):
        return hash(self.bond_dirs)

    def __repr__(self):
        return "Conformation(\"{}\")".format("".join(self.bond_dirs))

    def __str__(self):
        return "".join(self.bond_dirs)

    @property
    def locations(self):
        locations = [(self.location_delta[0], self.location_delta[1])]
        for bond_dir in self.bond_dirs:
            locations.append(next_monomer_location(locations[-1], bond_dir))
        return locations

    @property
    def overlapping(self):
        return len(self.bond_dirs) != len(set(self.locations)) - 1

    def _forward_contacts(self):
        locations = self.locations
        index_to_contacts = {i: [] for i in range(len(locations))}
        for i, location in enumerate(locations):
            for bond_dir in ['U', 'R', 'D', 'L']:
                adjacent_location = next_monomer_location(location, bond_dir)
                if adjacent_location in locations:
                    j = locations.index(adjacent_location)
                    if j > i + 1:
                        index_to_contacts[i].append(j)
        return index_to_contacts

    @cached_property
    def contacts(self):
        pairs = []
        for aa1, aa_forward_contacts in self._forward_contacts().items():
            for aa2 in aa_forward_contacts:
                pairs.append((aa1, aa2))
        return pairs

    def rotated_clockwise(self, n=1):
        rotated_bond_dirs = self.bond_dirs
        for _ in range(n):
            rotated_bond_dirs = [bond_dir_rotated_clockwise(bond_dir) for bond_dir in rotated_bond_dirs]
        return Conformation("".join(rotated_bond_dirs))

    @property
    def in_standard_form(self):
        i = 0
        bond_dir = self.bond_dirs[0]
        if bond_dir == 'U':
            return False
        while bond_dir == 'U':
            i += 1
            if not self.bond_dirs[i] in ['U', 'R']:
                return False
        return True

    @property
    def width(self):
        return 1 + abs(self.bond_dirs.count('R') - self.bond_dirs.count('L'))

    @property
    def height(self):
        return 1 + abs(self.bond_dirs.count('U') - self.bond_dirs.count('D'))

    def overlaps_with(self, other):
        return not set(self.locations).isdisjoint(other.locations)

    def contacts_with(self, other):
        contacts = []
        for i, location in enumerate(self.locations):
            for bond_dir in ['U', 'R', 'D', 'L']:
                adjacent_location = next_monomer_location(location, bond_dir)
                if adjacent_location in other.locations:
                    contacts.append((i, other.locations.index(adjacent_location)))
        return contacts


def generate_full_conformation_space(L):
    """Return all self-avoiding walks of length 'L' with the first bond Up and the first non-Up bond Right."""

    conformations = []
    n = L - 2  # index of last bond in 'conformation'
    first_R = n  # index of the first 'R' in the conformation
    conformation = ['U'] * (n + 1)
    while first_R > 0:
        # See if the current conformation has overlap
        x = y = j = 0
        res_positions = {(x, y): j}  # keyed by coords, items are residue numbers
        res_coords = [(x, y)]  # 'res_coords[j]' is coords of residue 'j'
        for c in conformation:
            x, y = next_monomer_location((x, y), c)
            if (x, y) in res_positions:  # overlap
                # increment at the step that gave the problem
                for k in range(j + 1, n + 1):
                    conformation[k] = 'U'
                conformation[j] = bond_dir_rotated_clockwise(conformation[j])
                while conformation[j] == 'U':
                    j -= 1
                    conformation[j] = bond_dir_rotated_clockwise(conformation[j])
                if j == first_R and conformation[j] not in ['R', 'U']:
                    first_R -= 1
                    conformation[first_R] = 'R'
                    for k in range(j, n + 1):
                        conformation[k] = 'U'
                break
            j += 1
            res_positions[(x, y)] = j
            res_coords.append((x, y))
        else:  # loop finishes normally, this is a valid conformation
            # generate the next conformation
            conformations.append(Conformation("".join(conformation)))
            i = n
            conformation[i] = bond_dir_rotated_clockwise(conformation[i])
            while conformation[i] == 'U':
                i -= 1
                conformation[i] = bond_dir_rotated_clockwise(conformation[i])
            # make sure first non-'U' is 'R'
            if i == first_R and conformation[i] not in ['R', 'U']:
                first_R -= 1
                conformation[first_R] = 'R'
                for j in range(i, n + 1):
                    conformation[j] = 'U'

    return conformations


class Ensemble:
    def __init__(self, conformations=None):
        if conformations is None:
            conformations = []

        self.contact_sets_to_conformations = dict()
        self.contact_sets = []
        for conformation in conformations:
            self.add(conformation)

    def add(self, conformation):
        contacts = frozenset(conformation.contacts)
        try:
            self.contact_sets_to_conformations[contacts].append(conformation)  # TODO: what if conformation changes?
        except KeyError:
            self.contact_sets_to_conformations[contacts] = [conformation]
            self.contact_sets.append(contacts)

    def conformations_with_contact_set(self, contact_set):
        return self.contact_sets_to_conformations[contact_set]

    @cached_property
    def contact_set_multiplicities(self):
        return [len(self.conformations_with_contact_set(cs)) for cs in self.contact_sets]

    @cached_property
    def num_conformations(self):
        return sum(self.contact_set_multiplicities)

    @property
    def conformations(self):
        return chain(*self.contact_sets_to_conformations.values())


class Lattice:
    def __init__(self, L, interaction_energies=miyazawa_jernigan):
        self.L = L
        self.interaction_energies = interaction_energies

    def __eq__(self, other):
        return self.L == other.L and self.interaction_energies == other.interaction_energies

    def __hash__(self):
        return hash(tuple([self.L, frozenset(self.interaction_energies.items())]))

    @cached_property
    def ensemble(self):
        return Ensemble(generate_full_conformation_space(self.L))

    def sum_contact_energy(self, seq, contact_set, seq2=None):
        seq2 = seq if seq2 is None else seq2
        energy = 0
        for i, j in contact_set:
            aa1 = seq[i]
            aa2 = seq2[j]
            energy += self.interaction_energies[aa1 + aa2]
        return energy  # TODO: precision handling

    def energy(self, seq, conformation):
        return self.sum_contact_energy(seq, conformation.contacts)

    @lru_cache(maxsize=16)
    def contact_set_energies(self, seq):
        n_contact_sets = len(self.ensemble.contact_sets)
        contact_set_energies = np.empty(n_contact_sets)
        for i, contact_set in enumerate(self.ensemble.contact_sets):
            contact_set_energies[i] = self.sum_contact_energy(seq, contact_set)
        return contact_set_energies

    def conformation_energies(self, seq):
        return np.repeat(self.contact_set_energies(seq), self.ensemble.contact_set_multiplicities)

    @lru_cache(maxsize=16)
    def minE_conformations(self, seq):
        energies = self.contact_set_energies(seq)
        minE_indices = np.where(energies == energies.min())[0]
        minE_contact_sets = [self.ensemble.contact_sets[i] for i in minE_indices]
        minE_conformations = []
        for contact_set in minE_contact_sets:
            for conformation in self.ensemble.conformations_with_contact_set(contact_set):
                minE_conformations.append(conformation)
        return minE_conformations

    def fold(self, protein):
        protein.conformations = self.minE_conformations(protein.seq)
        protein.lattice = self
        protein.energy = self.energy(protein.seq, protein.conformations[0])


class ProteinNotFoldedError(Exception):
    """Protein hasn't been folded yet, so certain of its attributes are not known."""


@dataclass
class Protein:
    seq: str
    conformations: List[Conformation] = None
    lattice: Lattice = None
    energy: float = None

    @property
    def native_state(self):
        if self.conformations is None:
            raise ProteinNotFoldedError

        return self.conformations[0] if len(self.conformations) == 1 else None

    @native_state.setter
    def native_state(self, conformation):
        self.conformations = [conformation]

    def partition_factor(self, temp=1.0):
        normalized_energies = np.append(-self.lattice.contact_set_energies(self.seq)/temp, -self.energy/temp)
        multiplicities = self.lattice.ensemble.contact_set_multiplicities + [-len(self.conformations)]
        return logsumexp(normalized_energies, b=multiplicities)

    def stability(self, temp=1.0):
        if self.conformations is None:
            raise ProteinNotFoldedError
        assert isinstance(temp, (int, float)) and temp > 0

        return self.energy + temp * self.partition_factor(temp)

    def frac_folded(self, temp=1.0):
        assert isinstance(temp, (int, float)) and temp > 0
        return 1.0 / (1.0 + np.exp(self.stability(temp) / temp))

    def bind(self, other):
        assert self.lattice.interaction_energies == other.lattice.interaction_energies
        ligand = deepcopy(other)
        min_binding_energy = inf
        min_binding_energy_rotations = []
        min_binding_energy_location_deltas = []

        unrotated_ligand_native_state = ligand.native_state
        for rotations in range(4):
            ligand.native_state = unrotated_ligand_native_state.rotated_clockwise(rotations)
            ligand_width = ligand.native_state.width
            ligand_height = ligand.native_state.height
            delta_x_range = range(-ligand_width, self.native_state.width + ligand_width)
            delta_y_range = range(-ligand_height, self.native_state.height + ligand_height)  # TODO: check this
            for delta_x in delta_x_range:
                for delta_y in delta_y_range:
                    ligand.native_state.location_delta = (delta_x, delta_y)
                    if self.native_state.overlaps_with(ligand.native_state):
                        continue
                    contact_set = self.native_state.contacts_with(ligand.native_state)
                    binding_energy = self.lattice.sum_contact_energy(self.seq, contact_set, ligand.seq)
                    if binding_energy < min_binding_energy:
                        min_binding_energy = binding_energy
                        min_binding_energy_rotations = [rotations]
                        min_binding_energy_location_deltas = [(delta_x, delta_y)]
                    elif binding_energy == min_binding_energy:
                        min_binding_energy_rotations.append(rotations)
                        min_binding_energy_location_deltas.append((delta_x, delta_y))

        return min_binding_energy, list(zip(min_binding_energy_rotations, min_binding_energy_location_deltas))
