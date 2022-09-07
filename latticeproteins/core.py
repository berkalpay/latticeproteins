import os
import math
import numpy as np
from copy import copy
from latticeproteins.interactions import miyazawa_jernigan


def next_monomer_location(location, bond_dir):
    bond_dir_to_dx = {'U': 0, 'R': 1, 'D': 0, 'L': -1}
    bond_dir_to_dy = {'U': 1, 'R': 0, 'D': -1, 'L': 0}
    return location[0] + bond_dir_to_dx[bond_dir], location[1] + bond_dir_to_dy[bond_dir]


class Conformation:
    def __init__(self, bond_dirs):
        self.bond_dirs = bond_dirs
        self._set_locations()

    def __getitem__(self, item):
        return self.bond_dirs[item]

    def __setitem__(self, key, value):
        self.bond_dirs[key] = value

    def __iter__(self):
        yield from self.bond_dirs

    def _set_locations(self):
        location = (0, 0)
        locations_to_index = dict()
        for i, bond_dir in enumerate(self.bond_dirs):
            locations_to_index[location] = i
            location = next_monomer_location(location, bond_dir)
        self.locations_to_index = locations_to_index

    def get_locations(self):
        return self.locations_to_index.keys()

    def overlapping(self):
        return len(self.bond_dirs) != len(set(self.get_locations()))

    def forward_contacts(self):
        index_to_contacts = dict()
        for i, location in enumerate(self.get_locations()):
            for bond_dir in ['U', 'R', 'D', 'L']:
                adjacent_location = next_monomer_location(location, bond_dir)
                if adjacent_location in self.get_locations():
                    j = self.locations_to_index[adjacent_location]
                    if j > i + 1:
                        index_to_contacts[i].append(j)
        return index_to_contacts


class Lattice:
    def __init__(self, L, interaction_energies=miyazawa_jernigan, pickle_dir="database/"):
        self.L = L
        self.interaction_energies = interaction_energies

        # Generate conformations
        conformations = []
        dx = {'U': 0, 'R': 1, 'D': 0, 'L': -1}
        dy = {'U': 1, 'R': 0, 'D': -1, 'L': 0}
        next = {'U': 'R', 'R': 'D', 'D': 'L', 'L': 'U'}
        n = self.L - 2  # index of last bond in 'conformation'
        first_R = n  # index of the first 'R' in the conformation
        conformation = Conformation(['U'] * (n + 1))
        while True:
            # See if the current conformation has overlap
            x = y = j = 0
            res_positions = {(x, y): j}  # keyed by coords, items are residue numbers
            res_coords = [(x, y)]  # 'res_coords[j]' is coords of residue 'j'
            for c in conformation:
                x += dx[c]
                y += dy[c]
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
                conformations.append(copy(conformation))
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
            # see if we are done
            if first_R == 0:
                break

        self.conformations = conformations

    def fold(self, protein):
        # Accelerate by Cython
        pass


class Protein:
    def __init__(self, seq, conformation=None):
        self.seq = seq
        self.conformation = conformation

    def set_conformation(self, conformation):
        self.conformation = conformation

    # TODO: should this be here?
    def energy(self, conformation):
        contacts = conformation.forward_contacts()
        energy = sum([interactions[c] for c in contacts])
        return energy
