from latticeproteins.core import Conformation

def test_locations():
    true_locs = {(0, 0), (0, 1), (0, 2), (0, 3), (1, 3), (1, 4),
                 (2, 4), (2, 3), (2, 2), (2, 1), (1, 1), (0, 1)}
    assert set(Conformation("UUURURDDDLL").get_locations()) == true_locs

def test_overlapping():
    assert Conformation("UUURURDDDLL").overlapping()

def test_non_overlapping():
    assert not Conformation("UUURURUUR").overlapping()

def test_contacts():
    assert len(Conformation("URDRUULLLDDDRRDLL").contacts()) == 10
