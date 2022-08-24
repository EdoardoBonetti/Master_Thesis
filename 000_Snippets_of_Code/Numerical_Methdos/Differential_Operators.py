####################################################
# Differential Operators
####################################################

from ngsolve import *

def Grad(u):
    # calculate the gradient of u with u = gf_u or u = cf_u
    if u.dim == 3:
        return CF( (u[0].Diff(x), u[1].Diff(y), u[2].Diff(z)) )
    if u.dim == 9:
        return CF( (Grad(u[0,:]),Grad(u[1,:]),Grad(u[2,:])),dims=(3,3) )

def Curl(u):
    # calculate the curl of u with u = gf_u or u = cf_u
    if u.dim == 3:
        return CF( (u[2].Diff(y)- u[1].Diff(z), u[0].Diff(z)- u[2].Diff(x), u[1].Diff(x)- u[0].Diff(y)) )
    if u.dim == 9:
        return CF( (Curl(u[0,:]),Curl(u[1,:]),Curl(u[2,:])),dims=(3,3) )   

def Div(u):
    # calculate the divergence of u with u = gf_u or u = cf_u
    if u.dim == 3:
        return CF( (u[0].Diff(x)+u[1].Diff(y)+u[2].Diff(z)) )
    if u.dim == 9:
        return CF( (Div(u[0,:]),Div(u[1,:]),Div(u[2,:])) )
    
def L2Norm(u, mesh):
    return sqrt(Integrate(InnerProduct(u,u), mesh))

def H1SemiNorm(u, mesh):
    return sqrt(Integrate(InnerProduct(Grad(u),Grad(u)), mesh))

def H1Norm(u, mesh):
    return sqrt(Integrate(InnerProduct(Grad(u),Grad(u)), mesh)) + L2Norm(u, mesh)

def HcurlNorm(u, mesh):
    return L2Norm(Curl(u), mesh) + L2Norm(u, mesh)

def HdivNorm(u, mesh):
    return L2Norm(Div(u), mesh) + L2Norm(u, mesh)

