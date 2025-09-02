import pysegtree
from pysegtree import *

def test_basic():
    print(SegmentTree.__dict__)
    segtree = SegmentTree([1.0, 2.0, 3.0])
    assert segtree is not None