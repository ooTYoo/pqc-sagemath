from sage.all import *
from mont_curve import *

def isogeny2(kerpt):
    gf = kerpt.gf
    tmp = kerpt.double()
    assert tmp.is_inf() == True
    #A' = 2(1-2x^2), B'=Bx
    A = gf(2)*(gf(1)-gf(2)*kerpt.x**2)
    B = kerpt.ec.b * kerpt.x
    ec = EC(A,B) # ec/<kerpt>
    #map pt of ec to ec/<kerpt>
    def map_of_pt(pt):
        if pt.x == kerpt.x:
            return EC_PT(kerpt.gf(0), kerpt.gf(0), ec)
        x = (kerpt.x * (pt.x **2) - pt.x)/(pt.x - kerpt.x)
        y = (pt.x**2)*kerpt.x - pt.gf(2)*pt.x*(kerpt.x**2) + kerpt.x
        y = y/(pt.x - kerpt.x)**2
        y = pt.y * y
        return EC_PT(x,y,ec)
    return (ec, map_of_pt)


def isogeny4(kerpt):
    gf = kerpt.gf
    tmp = kerpt.double().double()
    assert tmp.is_inf() == True
    #A' = 4x^4-2, B'= -x(x^2+1)*B/2
    A = gf(4)*(kerpt.x**4)-gf(2)
    B = -kerpt.ec.b * kerpt.x*(kerpt.x**2+gf(1)) /gf(2)
    ec = EC(A,B) # ec/<kerpt>

    #map pt of ec to ec/<kerpt>
    def map_of_pt(pt):
        t1 = pt.x - kerpt.x
        t2 = kerpt.gf(2)*pt.x*kerpt.x - (kerpt.x**2) - kerpt.gf(1)
        if t1*t2 == kerpt.gf(0):
            return EC_PT(kerpt.gf(0), kerpt.gf(0), ec)
        x = -pt.x*(pt.x*kerpt.x - kerpt.gf(1))**2
        x = x*(pt.x*(kerpt.x**2)+ pt.x -kerpt.gf(2)*kerpt.x)
        x = x/(t1**2)/t2
        y = -kerpt.gf(2)*(kerpt.x**2)*(pt.x*kerpt.x - kerpt.gf(1))
        y = y - kerpt.gf(4)*(pt.x**3)*(kerpt.x**3 + kerpt.x)
        y = y + kerpt.gf(2)*(pt.x**2)*(kerpt.x**4 + kerpt.gf(5)*kerpt.x**2)
        y = y - kerpt.gf(4)*pt.x*(kerpt.x**3 + kerpt.x)
        y = y + kerpt.x**2 + kerpt.gf(1)
        y = y/t1**3/t2**2
        y = pt.y * y
        return EC_PT(x,y,ec)
    return (ec, map_of_pt)

def isogeny3(kerpt):
    gf = kerpt.gf
    tmp = kerpt.triple()
    assert tmp.is_inf() == True
    #A'=(Ax-6x^3+6)*x, B'= Bx^2
    A = (kerpt.ec.a*kerpt.x-gf(6)*kerpt.x**2+gf(6))*kerpt.x
    B = kerpt.ec.b * kerpt.x**2
    ec = EC(A,B) # ec/<kerpt>

    #map pt of ec to ec/<kerpt>
    def map_of_pt(pt):
        t1 = pt.x - kerpt.x
        if t1 == kerpt.gf(0):
            return EC_PT(kerpt.gf(0), kerpt.gf(0), ec)
        t2 = kerpt.x*pt.x - kerpt.gf(1)
        x = pt.x * (t2**2)/(t1**2)
        y = (pt.x**2)*kerpt.x - kerpt.gf(3)*pt.x*(kerpt.x**2) + pt.x + kerpt.x
        y = pt.y * t2 * y
        y = y/(t1**3)
        return EC_PT(x,y,ec)
    return (ec, map_of_pt)

