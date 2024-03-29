#general idea of MPKC
import pdb
gf = GF(256,repr='int')
set_random_seed(0)
def gen_invertable_matrix(dim):
    mat = None
    while True:
        mat = Matrix([[gf.random_element() for i in range(dim)] for j in range(dim)])
        if mat.det() != gf(0):
            break
    inv = mat.inverse()
    return (mat,inv)

def gen_vect(dim):
    return vector([gf.random_element() for _ in range(dim)])

def gen_poly(n):
    coeff2 = Matrix(gf,n,n)
    coeff1 = vector([gf.random_element() for _ in range(n)])
    coeff0 = gf.random_element()
    for i in range(n):
        for j in range(i,n):
            coeff2[i,j] = gf.random_element()
    return (coeff2,coeff1,coeff0)

def gen_poly_ov(v,o,dim=None):    
    m = v + o
    if dim is None:
        dim = m
    coeff2 = Matrix(gf,dim,dim)
    coeff1 = vector(gf,dim)
    coeff0 = gf.random_element()
    for i in range(v):
        for j in range(i,v):
            coeff2[i,j] = gf.random_element()
        for k in range(o):
            coeff2[i,v+k] = gf.random_element()
    for i in range(m):        
            coeff1[i] = gf.random_element()
    return (coeff2,coeff1,coeff0)

def gen_cmap_mov(v,olist):
    #olist=[o1,o2], for rainbow there are 2 layers.
    assert(len(olist)==2)
    polys=[]
    n = v + sum(olist)
    for i in range(olist[0]):
        poly = gen_poly_ov(v,olist[0],dim=n)
        polys.append(poly)
    for i in range(olist[1]):
        poly = gen_poly_ov(v+olist[0],olist[1],dim=n)
        polys.append(poly)
    return polys

def upper_t_form(mat):
    m = mat.ncols()
    for i in range(m):
        for j in range(i+1,m):
            mat[i,j] += mat[j,i]
            mat[j,i] = gf(0)
    return mat

def keygen(v,olist):
    m = sum(olist)
    n = v + m
    t1 = gen_invertable_matrix(n)
    ct = gen_vect(n)
    s1 = gen_invertable_matrix(m)
    cs = gen_vect(m)
    cmap = gen_cmap_mov(v,olist)
    t,s = t1[0],s1[0]
    As,Bs,Cs = [],[],[]
    pub = []
    for poly in cmap:
        A,B,C = poly
        AA = t.transpose()*A*t
        AA = upper_t_form(AA)
        BB = B*t + ct*A*t + t.transpose()*A*ct
        CC = ct*A*ct + B*ct + C
        As.append(AA)
        Bs.append(BB)
        Cs.append(CC)
    for i in range(m):
        A,B,C = Matrix(gf,n,n),vector(gf,n),gf(0)
        for j in range(m):
            A += s[i,j]*As[j]
            B += s[i,j]*Bs[j]
            C += s[i,j]*Cs[j]
        pub.append((A,B,C + cs[i]))
    simiPK = [(As[i],Bs[i],Cs[i]) for i in range(m)]
    sk=(cmap,t1[1],ct,s1[1],cs)
    return (pub,sk,simiPK)

def poly_eval(poly,vect):
    A,B,C = poly
    return vect*A*vect + B*vect + C

def polys_eval(pk,vect):
    rslt = vector([poly_eval(p,vect) for p in pk])
    return rslt

def cmap2linear(polys,vect,tv):
    assert len(polys)==len(vect)
    o = len(vect)
    v = len(tv)
    m = o + v    
    coeff = Matrix(gf,o,o)
    bv = vector(gf,o)
    for i,poly in enumerate(polys):
        A,B,C = poly
        c0 = A[0:v,0:v]
        c1 = A[0:v,v:m]            
        coeff[i,:] = tv*c1 + B[v:m]
        bv[i] = vect[i] - tv*c0*tv - tv*B[0:v] - C
    return (coeff,bv)


def gauss_solv(AA,BB,tv):    
    if AA.det() == gf(0):
        return None
    rv = AA.inverse() * BB
    rslt = vector(list(tv) + list(rv))
    return rslt

def reverse_cmap(cmap,v,olist,vect):    
    o1,o2 = olist[0],olist[1]
    while True:
        tv1 = vector([gf.random_element() for _ in range(v)])
        AA,BB = cmap2linear(cmap[:o1],vect[:o1],tv1)
        ans1 = gauss_solv(AA,BB,tv1)
        if ans1 is None:
            continue
        AA,BB = cmap2linear(cmap[o1:],vect[o1:],ans1)
        ans2 = gauss_solv(AA,BB,ans1)
        if ans2 is None:
            continue
        return ans2

def rainbow_sign(h,sk, v, olist):
    F, invt,ct, invs, cs = sk    
    tmp = invs * (h - cs)
    #print("sign:",h)
    #print("sign: affin_inv_S:",tmp)    
    tmp = reverse_cmap(F,v,olist,tmp)
    #print("sign: affin_inv_F:",tmp)    
    sign = invt * (tmp - ct)
    #print("sign: affin_inv_T:",sign)    
    return sign

def verify_by_step(sk,psk,sign,h):
    F, invt,ct, invs, cs = sk
    #pdb.set_trace()
    t1 = invt.inverse() * sign + ct
    t2 = polys_eval(F,t1)
    t3 = invs.inverse() * t2 + cs
    print("verify input:",sign)
    print("verify affine_T:",t1)
    print("verify F:",t2)
    print("verify affine_S:",t3)
    tmp = polys_eval(psk,sign)
    print("test F:",tmp)
    calc = polys_eval(pk,sign)
    print("test h:",calc)
    return calc==h

def verify(pk,sign,h):
    calc = polys_eval(pk,sign)
    return calc==h


if __name__ == "__main__":
    v = 30
    olist = [15,15]
    pk, sk, hpk = keygen(v,olist)
    h = gen_vect(sum(olist))
    s = rainbow_sign(h,sk,v,olist)
    #verify_by_step(sk,hpk,s,h)
    r = verify(pk,s,h)
    print(r)
