class EC():
    def __init__(self,a,b):
        assert a.parent() == b.parent() #on the same fp
        self.a = a
        self.b = b
        self.gf = self.a.parent()   #the finit field
        jinv = self.gf(256)*(self.a**2- self.gf(3))**3/(self.a**2- self.gf(4))
        self.jinv = jinv

class EC_PT():
    def __init__(self,x,y,domain):
        assert x.parent()==y.parent()
        assert x.parent()==domain.gf
        self.x = x
        self.y = y
        self.ec = domain
        self.gf = self.ec.gf
        self.inf = False
        if self.x == self.gf(0) and self.y== self.gf(0):
            self.inf = True

    def is_inf(self):
        return self.inf

    def is_oncurve(self):
        if self.is_inf():
            return True
        # test B*y^2 = x^3 + A x^2 + x
        return self.ec.b*(self.y**2) ==  (self.x**3) + self.ec.a*(self.x**2) + self.x 

    def is_same_pt(self,pB):
        if self.ec.jinv == pB.ec.jinv:
            return self.x==pB.x and self.y==pB.y
        return False

    def deepcopy(self):
        return EC_PT(self.x,self.y,self.ec)

    def minus(self):
        return EC_PT(self.x, -self.y, self.ec)

    def double(self):
        if self.is_inf():
            return self.deepcopy()

        t1 = self.x**2 - self.gf(1)     #x^2-1
        t2 = (self.x**2 + self.ec.a*self.x + self.gf(1))*self.x     #x(x^2+Ax+1)
        if t2==self.gf(0):
            return EC_PT(x=self.gf(0), y=self.gf(0), domain=self.ec)

        xp1 = t1**2
        xp2 = self.gf(4)*t2
        x = xp1 / xp2
        # x^4 + 2Ax^3+6x^2+2Ax+1
        yp1 = self.x**4 + self.gf(2)*self.ec.a*(self.x**3) + self.gf(6)*(self.x**2) + self.gf(2)*self.ec.a*self.x + self.gf(1)
        yp1 *= t1
        yp2 = self.gf(8)* (t2**2)
        y = self.y * yp1 / yp2

        return EC_PT(x, y, domain=self.ec)


    def add(self, pB):
        assert self.ec == pB.ec
        if self.is_same_pt(pB):
            return self.double()
        if pB.is_inf():
            return self.deepcopy()
        if self.is_inf():
            return pB.deepcopy()
        if self.x == pB.x:
            return EC_PT(x=self.gf(0), y=self.gf(0),domain=self.ec)

        lmd = (self.y - pB.y) / (self.x - pB.x)
        #Br^2-(x+x_Q)-A
        x = self.ec.b* (lmd**2) - (self.x + pB.x) - self.ec.a
        y = lmd * (self.x - x) - self.y
        return EC_PT(x, y, domain=self.ec)

    def triple(self):
        if self.is_inf():
            return self.deepcopy()
        
        # x^4 - 4Ax -6x^2 -3
        t1 = (self.x**4) - self.gf(4)*self.ec.a*self.x - self.gf(6)*(self.x**2) -self.gf(3)
        # 4Ax^3 + 3x^4 + 6x^2 - 1
        t2 = self.gf(4)*self.ec.a*(self.x**3) + self.gf(3)*(self.x**4) + self.gf(6)*(self.x**2) - self.gf(1)
        if t2==self.gf(0):
            return EC_PT(x=self.gf(0), y=self.gf(0), domain=self.ec)

        x = t1**2 * self.x / (t2**2)
        t3 = (self.x**8) + self.gf(4)*self.ec.a*(self.x**7) + self.gf(28)*(self.x**6) + self.gf(28)*self.ec.a*(self.x**5)
        t3 += (self.gf(16)*(self.ec.a**2)+self.gf(6))*(self.x**4) + self.gf(28)*self.ec.a*(self.x**3) + self.gf(28)*(self.x**2)
        t3 += self.gf(4)*self.ec.a*self.x + self.gf(1)
        y = self.y * t1 * t3/ (t2**3)
        return EC_PT(x, y, domain=self.ec)

    def scalar_product(self, k):
        rslt = self.deepcopy()
        for i in bin(k)[3:]:
            rslt = rslt.double()
            if i=='1':
                rslt = self.add(rslt)
        return rslt