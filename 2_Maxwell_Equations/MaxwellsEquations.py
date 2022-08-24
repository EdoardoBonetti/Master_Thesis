from ngsolve import *
from netgen import gui
from netgen.csg import unit_cube
import matplotlib.pyplot as plt

def MaxwellPropagation(mesh , order = 1 , initail_E_field = None , bonus_intorder = 10 ,draw_e_field = True , draw_d_field = False, draw_divd_field = False , final_time  = 0.1, step_time = 0.002 , MoreInfo = False):
    
    # define the 3 spaces we need
    fes =  HCurl(mesh, order=order )*HDiv(mesh, order=order   , dirichlet= ".*")
    fes_projection =  HDiv(mesh, order=order   , dirichlet= ".*")*L2(mesh,  order=order-1)
    if MoreInfo:
        pass
    #electric field
    elfield = GridFunction(fes) 
    e, d = elfield.components
    e.Set ( initail_E_field, bonus_intorder=10, dual=True)

    # projectional part
    Pelfield = GridFunction(fes_projection) 
    Pd, Pq = Pelfield.components

    # trial and test function for H(curl) X H(div)
    uc, ud , = fes.TrialFunction()
    vc, vd , = fes.TestFunction()

    # create the inverse of the submatrix M_e 
    uc,vc = HCurl(mesh, order=order ).TnT()
    edofs = fes.Range(0)
    emb_e = Embedding(fes.ndof, edofs)
    inv_e = BilinearForm(uc*vc*dx).Assemble().mat.Inverse(inverse="sparsecholesky")

    #create the inverse of M_d 
    ud,vd = HDiv(mesh, order=order   , dirichlet= ".*").TnT()
    ddofs = fes.Range(1)
    emb_d = Embedding(fes.ndof, ddofs)
    inv_d = BilinearForm(ud*vd*dx).Assemble().mat.Inverse(inverse="sparsecholesky")

    # create the passage matrix M_de : H(curl) -> H(div)
    c = HCurl(mesh, order=order ).TrialFunction()
    v = HDiv(mesh, order=order   , dirichlet= ".*").TestFunction()
    M_de = BilinearForm(trialspace=HCurl(mesh, order=order ), testspace=HDiv(mesh, order=order   , dirichlet= ".*") , geom_free = True)
    M_de +=  InnerProduct(c, v)*dx
    M_de.Assemble()

    # create the passage matrix M_ed : H(curl) -> H(div)
    c = HCurl(mesh, order=order ).TestFunction()
    v = HDiv(mesh, order=order   , dirichlet= ".*").TrialFunction()
    M_ed = BilinearForm(testspace=HCurl(mesh, order=order ), trialspace=HDiv(mesh, order=order   , dirichlet= ".*"), geom_free = True)
    M_ed +=  InnerProduct(c, v)*dx
    M_ed.Assemble()

    # creating the matrix for the operator curl-curl:
    u = HCurl(mesh, order=order ).TestFunction()
    v = HCurl(mesh, order=order ).TrialFunction()
    K = BilinearForm(HCurl(mesh, order=order ))
    K += -InnerProduct( curl(u), curl(v) )*dx
    K.Assemble()

    # compose the 2 final matrices
    Mup = emb_e @ inv_e @ M_ed.mat @ emb_d.T
    Mdown = emb_d @ inv_d @ M_de.mat @ inv_e  @ K.mat @ emb_e.T 

    M =  Mdown  

    # Create the projection operator for the divergence free part


    u, p  = fes_projection.TrialFunction()
    v, q  = fes_projection.TestFunction()

    P = BilinearForm(fes_projection)
    P +=  v*u * dx
    P +=  div( v)*p * dx
    P +=  div( u)*q * dx
    P += -1e-10*p*q * dx

    F = BilinearForm(fes_projection)
    F +=  v*u * dx

    P.Assemble()
    F.Assemble()

    invP = P.mat.Inverse(inverse="sparsecholesky")

    if draw_e_field :
        scene_e = Draw(e, mesh,clipping=(0,0,1))
    if draw_divd_field :
        scene_divd = Draw(div(d), mesh,clipping=(0,0,1) ,draw_surf=False)
    if draw_d_field :
        scene_d = Draw(d, mesh,clipping=(0,0,1))
    
    t =   0.1
    tau = 0.002
    h = 0
    Energy = []
    Redraw(True)
    with TaskManager(pajetrace=10**8):
      while h < t :
        h+=tau

        #we use only one matrix of the form described in the above markdowm
        elfield.vec.data += tau * Mdown *elfield.vec
        elfield.vec.data += tau * Mup *elfield.vec

        #elfield.vec.data += tau* (inv @ Mm ) @ inv * elfield.vec
        #scene_d.Redraw()
        # input("stop")
        Pd.vec.data = d.vec
        Pelfield.vec.data =  invP*(F.mat*Pelfield.vec)
        d.vec.data = Pd.vec
        Energy.append(Integrate(InnerProduct(e,e),mesh))
        if draw_e_field :
            scene_e.Redraw()
        #print(h)
        if draw_d_field :
            scene_divd.Redraw()
        if draw_d_field :
            scene_d.Redraw()
    
    



final_time = 0.3
step_time = 0.002

Maxh = 0.13
order = 3
mesh = Mesh(unit_cube.GenerateMesh(maxh=Maxh))

# electric permettivity and magnetic permeability
eps = CF((1), dims=(1,1))
mu =  CF((1), dims=(1,1))

# define initial conditions for the wave
peak = exp(-35*( (x-0.5)**2 + (y-0.5)**2 +(z-0.5)**2))
initail_E_field = CF ( (0, 0, peak ), dims=(3,1))


MaxwellPropagation(mesh , order = order ,initail_E_field = initail_E_field , final_time  = final_time, step_time = step_time)