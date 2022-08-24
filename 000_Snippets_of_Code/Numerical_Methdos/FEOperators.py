from ngsolve import *
import time
from ngsolve import *
from ngsolve.webgui import Draw
from netgen.csg import unit_cube
from ngsolve.krylovspace import CGSolver


def Curl(u):
    if u.dim == 3:
        return CF( (u[1].Diff(z)- u[2].Diff(y), u[2].Diff(x)- u[0].Diff(z), u[0].Diff(y)- u[1].Diff(x)) )
    if u.dim == 9:
        return CF( (Curl(u[0,:]),Curl(u[1,:]),Curl(u[2,:])),dims=(3,3) )
def Inc(u):
    return Curl((Curl(u)).trans)

def IncOperator(mesh,order = 2, g_field = None,BND_Inc_g_field = None,bonus_intorder=10):

    fesHCurlCurl = HCurlCurl(mesh, order=order, dirichlet= ".*")

    gf_g = GridFunction(fesHCurlCurl) 
    gf_g.Set ( g_field, bonus_intorder=bonus_intorder, dual=True)

    gf_Inc_G = GridFunction(fesHCurlCurl)
    gf_Inc_G.Set(BND_Inc_g_field, BND,  bonus_intorder=bonus_intorder, dual=True)

    u,v = fesHCurlCurl.TnT()
    n = specialcf.normal(3)
    Q_n = Id(3) - OuterProduct(n,n) 
    n_cross_v = CF( (Cross(n,v[0,:]),Cross(n,v[1,:]),Cross(n,v[2,:])), dims=(3,3) )
    t1 =specialcf.EdgeFaceTangentialVectors(3)[:,0]
    t2 =specialcf.EdgeFaceTangentialVectors(3)[:,1]
    e = specialcf.tangential(3,True)
    n1 = Cross( t1, e)
    n2 = Cross( t2, e) 

    # Mass matrix
    a = BilinearForm(fesHCurlCurl)
    a += InnerProduct(u,v)*dx 

    # linear form induced by the metric gfG
    f = LinearForm(fesHCurlCurl)

    # thetrahedron inc part
    f += InnerProduct(gf_g.Operator("inc"),v)*dx        
    # faces part:
    f += ( InnerProduct(Q_n*n_cross_v, curl(gf_g).trans) + Cross(gf_g*n,n)*(curl(v)*n) )*dx(element_boundary=True)
    # Edges components: t'*v*C_n*n
    f += (gf_g[n1,e]*v[e,t1] - gf_g[n2,e]*v[e,t2])*dx(element_vb=BBND)
    
    a.Assemble()
    f.Assemble()

    r = f.vec.CreateVector()
    r.data = f.vec - a.mat * gf_Inc_G.vec
        #inverse = CGSolver(a.mat, pre.mat , printrates='\r', maxiter=500,tol=1e-9)
    gf_Inc_G.vec.data += a.mat.Inverse(freedofs=fesHCurlCurl.FreeDofs(),inverse="sparsecholesky") * r
        #gfCurlCurlG.vec.data += inver
    Draw(mesh, gf_Inc_G)
    return gf_Inc_G


mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))
peak = exp(-25*( (x-0.5)**2 + (y-0.5)**2 +(z-0.5)**2))
PEAK = CF ( (peak, 0, 0 , 0,peak,0,0,0,peak), dims=(3,3))
IncPeak = CF( Inc(PEAK), dims=(3,3) )
bonus_intorder = 10
gf_Inc_G = IncOperator(mesh,order = 2, g_field =PEAK, BND_Inc_g_field = IncPeak ,bonus_intorder = 10)

Draw(mesh, gf_Inc_G)

input("nin")