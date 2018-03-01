"""This package contains utilities for lattice protein conformations.

Uses 2-dimensional non-compact models.

Originally written by Jesse Bloom.

Extended by Zach Sailer."""
from .sequences import random_sequence
from .protein import LatticeProtein, LatticeProteins
from .conformations import Conformations, ConformationList
from .interactions import miyazawa_jernigan
from . import draw
