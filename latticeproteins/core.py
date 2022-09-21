import math
from functools import cached_property
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
    dx = {'U': 0, 'R': 1, 'D': 0, 'L': -1}
    dy = {'U': 1, 'R': 0, 'D': -1, 'L': 0}
    next = {'U': 'R', 'R': 'D', 'D': 'L', 'L': 'U'}
    n = L - 2  # index of last bond in 'conformation'
    first_R = n  # index of the first 'R' in the conformation
    conformation = ['U'] * (n + 1)
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
        # see if we are done
        if first_R == 0:
            break

    return conformations


class Ensemble:
    def __init__(self, conformations=None):
        self.contact_sets_to_conformations = dict()
        if conformations:
            for conformation in conformations:
                self.add(conformation)

    def __iter__(self):
        yield from self.contact_sets_to_conformations

    def __getitem__(self, item):
        return self.contact_sets_to_conformations[item]

    def add(self, conformation):
        contacts = frozenset(conformation.contacts)
        try:
            self.contact_sets_to_conformations[contacts].append(conformation) # TODO: what if conformation changes?
        except KeyError:
            self.contact_sets_to_conformations[contacts] = [conformation]


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
        return energy

    def energy(self, seq, conformation):
        return self.sum_contact_energy(seq, conformation.contacts)

    def energies(self, seq):
        conformation_to_energy = dict()
        for contact_set in self.ensemble:
            energy = self.sum_contact_energy(seq, contact_set)
            for conformation in self.ensemble[contact_set]:
                conformation_to_energy[conformation] = energy
        return conformation_to_energy

    def minE_conformations(self, seq):
        minE = math.inf
        minE_conformations = []
        for conformation, energy in self.energies(seq).items():
            if energy < minE:
                minE = energy
                minE_conformations = [conformation]
            elif energy == minE:
                minE_conformations.append(conformation)
        return minE_conformations

    def fold(self, protein):
        minE_conformations = self.minE_conformations(protein.seq)
        if len(minE_conformations) == 1:
            protein.set_conformation(minE_conformations[0])


class Protein:
    def __init__(self, seq, conformation=None):
        self.seq = seq
        self.conformation = conformation

    def __repr__(self):
        return "Protein(seq={}, conformation={})".format(self.seq, self.conformation)

    def set_conformation(self, conformation):
        self.conformation = conformation
