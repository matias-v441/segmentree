import pysegtree
from pysegtree import *

def test_basic():
    print(SegmentTree.__dict__)
    segtree = SegmentTree([1.0, 2.0, 3.0])
    assert segtree is not None
    segtree.add_segment((1.0, 2.0), 0)
    segtree.add_segment((2.0, 3.0), 1)
    union = segtree.get_union((float('-inf'), float('inf')))
    assert union.contains_point(1.5)
    assert union.contains_point(2.5)
    assert not union.contains_point(3.5)
    s = segtree.root_stats
    assert s["length"] == 2.0
    assert s["max_ovp"] == 1
    assert s["min_ovp"] == 0