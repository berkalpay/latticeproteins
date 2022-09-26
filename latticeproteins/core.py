from functools import cached_property, lru_cache
from dataclasses import dataclass
from typing import List
from itertools import chain

import numpy as np
from scipy.special import logsumexp

from latticeproteins.interactions import miyazawa_jernigan


def next_monomer_location(location, bond_dir):
    bond_dir_to_dx = {'U': 0, 'R': 1, 'D': 0, 'L': -1}
    bond_dir_to_dy = {'U': 1, 'R': 0, 'D': -1, 'L': 0}
    return location[0] + bond_dir_to_dx[bond_dir], location[1] + bond_dir_to_dy[bond_dir]


class Conformation:
    def __init__(self, bond_dirs):
        self._bond_dirs = bond_dirs

    @property
    def bond_dirs(self):
        return self._bond_dirs

    @bond_dirs.setter
    def bond_dirs(self, value):
        self._bond_dirs = value
        self.__dict__.pop("contacts", None)

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

    def get_bond_dirs(self):
        return self.bond_dirs

    def set_bond_dirs(self, bond_dirs):
        self.bond_dirs = bond_dirs

    def get_locations(self):
        locations = [(0, 0)]
        for bond_dir in self.bond_dirs:
            locations.append(next_monomer_location(locations[-1], bond_dir))
        return locations

    def overlapping(self):
        return len(self.bond_dirs) != len(set(self.get_locations())) - 1

    def _forward_contacts(self):
        locations = self.get_locations()
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


def generate_full_conformation_space(L):
    # Generate conformations
    conformations = []
    next = {'U': 'R', 'R': 'D', 'D': 'L', 'L': 'U'}
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
                conformation[j] = next[conformation[j]]
                while conformation[j] == 'U':
                    j -= 1
                    conformation[j] = next[conformation[j]]
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
            conformation[i] = next[conformation[i]]
            while conformation[i] == 'U':
                i -= 1
                conformation[i] = next[conformation[i]]
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

    def conformations(self):
        return chain(*self.contact_sets_to_conformations.values())

class Lattice:
    def __init__(self, L, interaction_energies=miyazawa_jernigan):
        self.L = L
        self.interaction_energies = interaction_energies

    @cached_property
    def ensemble(self):
        return Ensemble(generate_full_conformation_space(self.L))

    def sum_contact_energy(self, seq, contact_set):
        energy = 0
        for i, j in contact_set:
            aa1 = seq[i]
            aa2 = seq[j]
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

    def partition_factor(self, temp=1.0):
        minE = self.lattice.energy(self.seq, self.conformations[0])
        conformation_energies = self.lattice.conformation_energies(self.seq)
        # TODO: optimize?
        return logsumexp(np.append(-self.lattice.conformation_energies(self.seq)/temp, -minE/temp * len(self.conformations)),
                         b=[1]*len(conformation_energies) + [-1])

    def stability(self, temp=1.0):
        if self.conformations is None:
            raise ProteinNotFoldedError
        assert isinstance(temp, (int, float)) and temp > 0

        return self.energy + temp * self.partition_factor(temp)

    def frac_folded(self, temp=1.0):
        assert isinstance(temp, (int, float)) and temp > 0
        return 1.0 / (1.0 + np.exp(self.stability(temp) / temp))
