#!/usr/bin/env sage
import sys
from sage.all import *
import sike_para
import mont_curve
import time
import isogeny

preparser(True)
para = sike_para.SIKEp434
lA = int(para["lA"], 10)
eA = int(para["eA"], 10)
lB = int(para["lB"], 10)
eB = int(para["eB"], 10)
a = int(para["A"])
b = int(para["B"])
p = int(para["p"],16)
gf = GF(p)
gfp.<x> = gf[]
gf2 = gf.extension(x^2+1,'i')

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

def gf2_pf(a):
    abar = a.conjugate()
    im_a = (a-abar)/2
    print("Re:",hex(int(a.polynomial()(0))))
    print("Im:", hex(int(im_a.polynomial()(1))))

def isogenkey(sk,party):
    '''
    sk : secret parameters to calc kerpt = Pl + [sk]Ql
    party: 'A' for 2, 'B' for 3
    output: (iosP,iosQ), which contians final curve=EC/ker,P and Q are basepoint of peer party
    '''
    if party=='A':
        P,Q,e,l,isomap = PA,QA,eA,lA,isogeny.isogeny2
    elif party=='B':
        P,Q,e,l,isomap = PB,QB,eB,lB,isogeny.isogeny3
    else:
        print("[!] wrong party-para, check input again!")
        return None
    S = Q.scalar_product(sk)
    S = S.add(P)
    # should map the basepoint of peer party to their isogeny image
    if party=='B':
        P,Q = PA,QA
    elif party=='A':
        P,Q = PB,QB

    for i in range(e):
        ord = l**(e-1-i)
        kerpt = S.scalar_product(ord)
        ec,ptmap = isomap(kerpt)
        P = ptmap(P)
        Q = ptmap(Q)
        S = ptmap(S)
    return (P,Q)

def isoex(sk,party,oppside_pk):
    '''
    sk : secret parameters to calc kerpt = opP + [sk]opQ
    party: 'A' for 2, 'B' for 3
    output: final cureve (EC), which contains EC.jinv
    '''
    if party == 'A':
        e, l, isomap = eA, lA, isogeny.isogeny2
    elif party == 'B':
        e, l, isomap = eB, lB, isogeny.isogeny3
    else:
        print("[!] wrong party-para, check input again!")
        return None
    P,Q = oppside_pk
    S = Q.scalar_product(sk)
    S = S.add(P)
    for i in range(e):
        ord = l ** (e - 1 - i)
        kerpt = S.scalar_product(ord)
        ec, ptmap = isomap(kerpt)
        S = ptmap(S)
    return S.ec

def test_para():
    #prime number test
    pp = (lA**eA)*(lB**eB)-1
    assert pp == p
    print("[+] test <P ==lA^eA * lB^eB - 1> --- OK")

    #on_cureve test
    assert QA.is_oncurve() == True
    assert PA.is_oncurve() == True
    assert QB.is_oncurve() == True
    assert QB.is_oncurve() == True
    print("[+] on cuver test of PA,QA,PB,QB --- OK")

    #ecc point calc test
    pt2 = QB.double()
    pt3 = QB.add(pt2)
    pt33 = QB.triple()
    pt333 = QB.scalar_product(3)
    assert pt3.is_same_pt(pt33)== True
    assert pt3.is_same_pt(pt333)==True
    print("[+] test <Point-double, Point-Add, Point-Triple, Point-Scalar-Product> --- OK")

    #para test
    ordA = lA**eA
    ordB = lB**eB
    assert PA.scalar_product(ordA).is_inf() == True
    assert QA.scalar_product(ordA).is_inf() == True
    assert PB.scalar_product(ordB).is_inf() == True
    assert QB.scalar_product(ordB).is_inf() == True
    print("[+] test order of base Point --- OK")

def sidh():
    beg = time.time()
    sk_A = int('caffbabedeadbeef',16)
    sk_B = int('caffbabedeadbeef'[::-1],16)
    pk_A = isogenkey(sk_A,'A')
    pk_B = isogenkey(sk_B,'B')
    ec1 = isoex(sk_A,'A',pk_B)
    ec2 = isoex(sk_B,'B',pk_A)
    print('keygen for A')
    #gf2_pf(ec1.a)
    #gf2_pf(ec1.b)
    gf2_pf(ec1.jinv)
    print('keygen for B')
    #gf2_pf(ec2.a)
    #gf2_pf(ec2.b)
    gf2_pf(ec2.jinv)
    end = time.time()
    print("key exchange -- ", 'OK' if ec1.jinv==ec2.jinv else 'Fail')
    print("[+] time consuming: ", end-beg, "s")


if __name__ == "__main__":
    test_para()
    sidh()