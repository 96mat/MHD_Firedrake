from firedrake import *
def ic_velocity_regular(x,y):
    u_ic=as_vector((y*exp(-4*(x**2+y**2)),\
               -x*exp(-4*(x**2+y**2))))
    return u_ic

def ic_magnetic_regular(x,y):
    B_ic_onZ=as_vector((2*x*y*sin(pi*y)-x*pi*(1-y**2)*cos(pi*y)+(1-x**2)*sin(pi*x),\
                    (1-y**2)*sin(pi*y)+2*x*y*sin(pi*x)-y*pi*(1-x**2)*cos(pi*x)))
    return B_ic_onZ


def ic_velocity_corona(x,y,xc,yc,r0,sigma):
    xc=Constant(xc)
    yc=Constant(yc)
    r0=Constant(r0)
    sigma=Constant(sigma)                  # Corona radius
    r  = sqrt((x - xc)**2 + (y - yc)**2)   # Thickness (smaller this value= thinner)
    g = exp(-((r - r0)**2) / (2*sigma**2)) # Anual Envelope
    u_ic_2 = g * as_vector((-(y - yc), (x - xc)))
    return u_ic_2


def ic_magnetic_Lshape(x,y):
    B_ic_onZ=as_vector((-x*(1-y**2)*sin(pi*y),\
                    -(1-x**2)*y*sin(pi*x)))
    return B_ic_onZ