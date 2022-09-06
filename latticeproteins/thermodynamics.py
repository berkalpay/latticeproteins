"""Module for calculating thermodynamics of lattice protein sequences."""

import os
import numpy as np

from . import conformations as c
from .interactions import miyazawa_jernigan
from .conformations import fold_energy


class LatticeThermodynamics(object):
    """A single lattice protein.

    Parameters
    ----------
    temp : float
        the temperature at which the fitness is computed.

    conformations : conformations.Conformations object
        is the 'conformations.Conformations' object
        used to fold the protein sequences.  'conformations.Length()'
        specifies the length of the protein sequences that can be
        folded.
    """

    def __init__(self, conformations, temp=1.0):
        assert isinstance(temp, (int, float)) and temp > 0
        assert isinstance(conformations, c.Conformations) or isinstance(conformations, c.ConformationList)

        self.temp = temp
        self.conformations = conformations

    def length(self):
        """Returns the sequence length for which fitnesses are computed."""
        return self.conformations.length()

    def _nativeE(self, seq):
        """Compute the lattice native energy and partition sum of a sequence.

        Parameters
        ----------
        seq : str or list
            sequence to fold.

        Returns
        -------
        minE : float
            Energy of the native state.
        conf : str
            Native conformation
        partitionsum : float
            Partition function sum
        numcontacts : int
            Number of contacts in native state
        folds : boolean
            True if a single native structure exists. False is not.
        """
        assert len(seq) == self.length()
        return self.conformations.fold_sequence(seq, self.temp)

    def stability(self, seq, target=None):
        """Computes the stability of a sequence if it is below cutoff.

        Parameters
        ----------
        seq : str or list
            sequence to fold.

        Returns
        -------
        stability : float
            Folding stability of the native state.
        """
        nativeE_results = self._nativeE(seq)
        if target is not None:
            minE = fold_energy(seq, target, self.conformations._interaction_energies)
            nativeE_results = list(nativeE_results)
            nativeE_results[0] = minE
            nativeE_results[1] = target
            nativeE_results[-1] = True
        stability_results = self._stability(*nativeE_results)
        return stability_results[0]

    def _stability(self, minE, conf, partitionsum, folds):
        """Computes a stability from minE and partition function.

        Parameters
        ----------
        seq : str or list
            sequence to fold.

        Returns
        -------
        stability : float
            Folding stability of the native state.
        conf : str
            Native conformation
        partitionsum : float
            Partition function sum
        numcontacts : int
            Number of contacts in native state
        folds : boolean
            True if a single native structure exists. False is not.
        """
        #(minE, conf, partitionsum, numcontacts, folds) = self._NativeE(seq)
        # Calculate a stability... if calculation does not work, stability = 0
        if folds:
            gu = - self.temp * np.log(partitionsum - np.exp(-minE / self.temp))
            dGf = minE - gu
            return (dGf, conf, partitionsum, folds)
        else:
            return (0, conf, partitionsum, folds)

    def fracfolded(self, seq, target=None):
        """Compute the fraction folded of the sequence.

        Parameters
        ----------
        seq : str or list
            sequence to fold.

        Returns
        -------
        fracfolded : float
            fractioned folded.
        """
        nativeE_results = self._nativeE(seq)
        if target is not None:
            minE = fold_energy(seq, target, self.conformations._interaction_energies)
            nativeE_results = list(nativeE_results)
            nativeE_results[0] = minE
            nativeE_results[1] = target
            nativeE_results[-1] = True
        stability_results = self._stability(*nativeE_results)
        fracfolded_results = self._fracfolded(*stability_results)
        return fracfolded_results[0]

    def _fracfolded(self, dG, conf, partitionsum, folds):
        """Computes the fitness from a given stability value

        Parameters
        ----------
        seq : str or list
            sequence to fold.

        Returns
        -------
        fracfolded : float
            Energy of the native state.
        conf : str
            Native conformation
        partitionsum : float
            Partition function sum
        numcontacts : int
            Number of contacts in native state
        folds : boolean
            True if a single native structure exists. False is not.
        """
        # folding to a target conformation
        #(dG, conf, partitionsum, numcontacts, folds) = self._Stability(seq)
        if folds is False:
            return (0, conf, partitionsum, folds)
        else:
            f = 1.0 / (1.0 + np.exp(dG / self.temp))
            return (f, conf, partitionsum, folds)
        

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
