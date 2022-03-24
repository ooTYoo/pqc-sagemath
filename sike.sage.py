#!/usr/bin/env sage
import sys
from sage.all import *
import sike_para
import mont_curve
import time

preparser(True)

para = sike_para.SIKEp434
p = int(para["p"],16)
gf = GF(p)
gfp.<x> = gf[]
gf2 = gf.extension(x^2+1,'i')
ec = None

a = int(para["A"])
b = int(para["B"])
ec = mont_curve.EC(gf2(a),gf2(b))
xQA = gf2(int(para["xQA0"], 16)+int(para["xQA1"], 16)*x)
yQA = gf2(int(para["yQA0"], 16)+int(para["yQA1"], 16)*x)
xPA = gf2(int(para["xPA0"], 16)+int(para["xPA1"], 16)*x)
yPA = gf2(int(para["yPA0"], 16)+int(para["yPA1"], 16)*x)
xQB = gf2(int(para["xQB0"], 16)+int(para["xQB1"], 16)*x)
yQB = gf2(int(para["yQB0"], 16)+int(para["yQB1"], 16)*x)
xPB = gf2(int(para["xPB0"], 16)+int(para["xPB1"], 16)*x)
yPB = gf2(int(para["yPB0"], 16)+int(para["yPB1"], 16)*x)
QA = mont_curve.EC_PT(xQA,yQA,ec)
PA = mont_curve.EC_PT(xPA,yPA,ec)
QB = mont_curve.EC_PT(xQB,yQB,ec)
PB = mont_curve.EC_PT(xPB,yPB,ec)

def test_para():
    global p,para
    la = int(para["lA"],10)
    ea = int(para["eA"],10)
    lb = int(para["lB"],10)
    eb = int(para["eB"],10)
    pp = (la**ea)*(lb**eb)-1
    assert pp == p
    print(hex(pp))

def gf2_pf(a):
    abar = a.conjugate()
    im_a = (a-abar)/2
    print("Re:",hex(int(a.polynomial()(0))))
    print("Im:", hex(int(im_a.polynomial()(1))))

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


if __name__ == "__main__":
    test_para()
    test_mont_curve()
