"""Module for calculating thermodynamics of lattice protein sequences."""

import numpy as np


def stability(E, partitionsum, temp):
    assert isinstance(temp, (int, float)) and temp > 0
    return E + temp * np.log(partitionsum - np.exp(-E / temp))


def frac_folded(stability, temp):
    assert isinstance(temp, (int, float)) and temp > 0
    return 1.0 / (1.0 + np.exp(stability / temp))
    # TODO: confirm the energies returned by ContactLooper the same as those returned by fold_energy?
