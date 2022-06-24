#general idea of MPKC

gf = GF(256)

def gen_invertable_matrix(dim):
	mat = None
	while ture:
		mat = [[gf.random_element() for i in range(dim)] for j in range(dim)]
		if mat.det() != gf(0):
			break
	inv = mat.inverse()
	return (mat,inv)

def gen_vect(dim):
	return [gf.random_element() for _ in range(dim)]

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
	coeff1 = vector([gf.random_element() for _ in range(m)])
	return (coeff2,coeff1,coeff0)

def gen_cmap_mov(v,olist):
	#olist=[o1,o2], for rainbow there are 2 layers.
	assert(len(olist)==2)
	polys=[]
	n = v + sum(olsit)
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
			C += s[i,j]*Cs[j] + cs[j]
		pub.append((A,B,C))

	sk=(cmap,t[1],ct,s[1],cs)
	return (pub,sk)






