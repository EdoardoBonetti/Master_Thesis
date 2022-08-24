from cmath import sqrt
from ngsolve import *

def InverseBlockDiagonal(fesU) :

    # from a Mass block matrix 
    #       |M_0    0 |
    #       | 0    M_1|
    # extract the inverse of matrix M_i for i in 0,1

    #fes = fesU*fesV
    u,v = fesU.TnT()
    
    print("OKKEY")
    return BilinearForm(u*v*dx).Assemble().mat.Inverse(inverse="sparsecholesky")


def MixedBlockMatrix(fesU,fesV) :
    
    # from a Mass block matrix 
    #       |M_00  M_01 |
    #       |M_10  M_11 |
    # 
    # extract the inverse of matrix M_i,(1-i) for i in 0,1 and j the 
    # i is the space of codomain 


    # creates the mixed matrix M_vu : fesU -> fesV
    u = fesV.TrialFunction()
    v = fesU.TestFunction()
    M_vu = BilinearForm(trialspace=fesU, testspace=fesV, geom_free = True)
    M_vu +=  InnerProduct(u, v)*dx
 
    return M_vu.Assemble()


