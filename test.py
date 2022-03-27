#!/usr/bin/env sage
import sys
from sage.all import *
import sike_para
import mont_curve
import time
from sike.sage.py import *

def test_para():
    pp = (lA**eA)*(lB**eB)-1
    assert pp == p
    print("[+] test <P ==lA^eA * lB^eB - 1> OK")

def test_mont_curve():
    assert QA.is_oncurve() == True
    assert PA.is_oncurve() == True
    assert QB.is_oncurve() == True
    assert QB.is_oncurve() == True
    pt2 = QA.double()
    pt3 = QA.add(pt2)
    pt33 = QA.triple()
    pt333 = QA.scalar_product(3)
    assert pt3.is_same_pt(pt33)== True
    assert pt3.is_same_pt(pt333)==True
    print("[+] test <Point-double, Point-Add, Point-Triple, Point-Scalar-Product> OK")

if __name__ == "__main__":
    test_para()
    test_mont_curve()
