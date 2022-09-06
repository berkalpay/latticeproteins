"""Module for calculating thermodynamics of lattice protein sequences."""

import numpy as np
import conformations as c


def stability(E, partitionsum, temp):
    assert isinstance(temp, (int, float)) and temp > 0
    return E + temp * np.log(partitionsum - np.exp(-E / temp))


def frac_folded(stability, temp):
    assert isinstance(temp, (int, float)) and temp > 0
    return 1.0 / (1.0 + np.exp(stability / temp))
    # TODO: confirm the energies returned by ContactLooper the same as those returned by fold_energy?


class LatticeGroupThermodynamics(object):
    """A group of lattice proteins. Useful for fast calculation of a bunch of
    lattice proteins with the same length.

    Parameters
    ----------
    seqlist : list
        List of lattice proteins.

    temp : float
        temperature of the system.

    conformations : Conformations or ConformationList object
        Conformation database for lattice with set length

    target : str (optional, default=None)
        target conformation to fold protein list

    Attributes
    ----------
    seqlist : list
        list of sequences.

    temp : float
        temperature of the system.

    nativeEs : array
        native (or target) energy for sequences in seqlist

    stabilities : array
        array of stabilities for sequences in seqlist

    fracfolded : array
        array of fraction folded for sequences in seqlist
    """
    def __init__(self, seqlist, temp, conformations, target=None):
        # Assign class instance variables and error check
        assert isinstance(temp, (int, float)) and temp > 0
        self.temp = temp

        # Set conformations after checking
        if not isinstance(conformations, c.Conformations) and not isinstance(conformations, c.ConformationList):
            raise TypeError("conformations must be a Conformations or ConformationList object.")
        self.conformations = conformations
        self.length = self.conformations.length

        # Check length of target
        assert len(target)-1 == self.length
        self._target = target

        self.seqlist = seqlist
        self.nativeEs = np.empty(len(self.seqlist), dtype=float)
        self.conformations = np.empty(len(self.seqlist), dtype="U|" + self.length)
        self._partitionsum = np.empty(len(self.seqlist), dtype=float)

        for i, seq in enumerate(self.seqlist):
            # Check that sequence is valid
            assert len(seq) == self.length

            # Fold sequence and set energy
            out = self.conformations.fold_sequence(seq, self.temp)
            if target is not None:
                self.nativeEs[i] = conformations.fold_energy(seq, target)
            else:
                self.nativeEs[i] = out[0]
            self.conformations[i] = out[1]
            self._partitionsum[i] = out[2]

    @property
    def stabilities(self):
        """Folding stability for all sequences in seqlist."""
        gu = - self.temp * np.log(self._partitionsum - np.exp(-self.nativeEs / self.temp))
        dGf = minE - gu
        return dGf

    @property
    def fracfolded(self):
        """Fracfolded folded for all sequences in seqlist."""
        return 1.0 / (1.0 + np.exp(self.stability / self.temp))
