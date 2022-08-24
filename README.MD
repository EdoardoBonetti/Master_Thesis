### AUTHORS :

    Edoardo Bonetti

### AKNOWLEDGMENT :

    Huge thanks to Michael Neunteufel and Joachim Schoeberl for helping me with both programming and theory behind this Prioject
    A thanks also to Matthias Hochsteger that gave me help with the code for creating the MP4 files

### IMPOROVEMNTS :

    There is the possibility to review the same work that I did here with the help of the new developed method created in NGSolve by J. Schoeberl and Markus (??) using the mixed-intersection triangle-honeycomb for the H(div)-H(curl), in the future the same work will be further divided to create such advances tecniques 

### DISCLAIMER :

    In thesource problem I divide between the direct problem and the inverse problem!

### IMPORTANT REMARKS and NOTATIONS :

1. All the functions are CoefficientFunctions are defined as `cf_name`.

2. If we apply a differential operator that we defined 
    ```python
        def DiffOpName(name):
            # name can be gf_name or cf_name
            ...
            return CoefficienFunction(...)
    ```
    we call the result `cf_DiffOpName_name = DiffOpName(cf_name)`.

3. All the gridfunctions are defined as `gf_name`.

4. If we apply a differential operator defined above then  we call the result `gf_DiffOpName_name = DiffOpName(gf_name)`.
    Even if the result is contained in a Coefficient Function, htis is done to reflect the possible discontinuity that the initial grid function may occur in if we keep differenciating.
    Also in the future it may be that we decide to change the definition such that if I enter a grid function, then it will return a grid function (same thing for a coefficient function to coefficient function). 

5. the the boundary conditions are defined as follows:
    
    i) in case of a dirichlet boundary condition:    

    `d_bc_cf_name` is a CoefficientFunction.

    `d_bc_gf_name` is a gridfunction.

    ii) the same holds for the Neumann boundary condition:
        
    `n_bc_cf_name` is a CoefficientFunction.
    
    `n_bc_gf_name` is a gridfunction.

    Obviously `name` could be `DiffOpName_name`. 

6. Throughout the code, when we want to invert a matrix, we use the comand
    `a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")` , notice that the flag `inverse = "pardiso"` can be substitute with `inverse = "sparsecholesky"` or, if NGSolve was compiled with the correct flags, we can also use `inverse = "(??)"`.

7. The general theory of mixed finite element is usually written as follows: :

    
    $$
    \text{Find the solution (u,p) of the problem:} \\
    \begin{cases}
    A(u,v) + B(v,p) = f(v) & \forall v\in V\\
    B(u,q) =g(q) & \forall q\in Q\\
    \end{cases}\\
    $$
    
    therefore we try to stick with the notation of test and trial functions even in the code:
    ```python
        (u,p) , (v,q) = fes.TnT() 
    ```
    but in specific cases where there is possibility of confusion then the best thing to do is to use the notation:
    ```python
        ts_name , tl_name  = fes.TnT() 
    ```
    and obviously the various subcomponents are called similarly

# Division of the Folder : `NGSOLVE_FOR_THESIS`

## `Books Thesis` :
    contain all the book for the thesis; If some book is needed in some code then we create the necessary "link" to it in the README of the folder

## `Curl_Problems` :
Study of the 3-dimensional curl problem with an intro of the 3x3-dimensional curl 

1. contains the `Helmotz_Decomposition.ipynb` of a function into the rotational and irrotational part (also the potential vectors are )
2. contains the `Source_Problem.ipynb` related to the curl operator
3. contains the `Initial_Problem.ipynb` related to the curl operator (TO DO) 
4. `Matrix_Curl_Problems` : 
subfolder containing a first on the Matrix Curl Problem in different FESpaces
    i. `Matrix_Curl_Problems.ipynb`


## Derivations :
1. `Derivation.ipynb` contains the derivations in Exact Math of some operators 

## Inc_Problems :
Very similar to the Curl Problem as structure

1. contains the `Beltrami_Helmotz_Decomposition.ipynb` of a function into the rotational and irrotational part (also the potential vectors are )
2. contains the `Source_Problem.ipynb` related to the Inc operator

3. contains the `Initial_Problem.ipynb` related to the Inc operator (TO DO) 
4. `Inc_Operator` folder containing :
    i. `2D_Inc_Operator.ipynb` a study of the 2 dimensional Inc operator.
        ii. `3D_Inc_Operator.ipynb` a study of the 3 dimensional Inc operator.

## `Einstein_Field_Equations` :

    We use the theoy of the `Inc_Problems` folder to Study the einstein equations

    1. `Source_Problem.ipynb`
    2. `Initial_Problem.ipynb`

## `Einstein_Field_Equations` :

    We use the theoy of the 'Curl_Problems' folder to Study the Maxwell's equations, in particular we only focus on the wave propagation odf the Electric wave equation:

    1. `Source_Problem.ipynb`
    2. `Initial_Problem.ipynb`
    3. `HCurl_L2_Projection.ipynb` : contains the projection needed to pass from ...


## `Snippets_Of_Code` :

    Contains all the snippets that might be usefull or that we use in the code explained more in details, for example how to extract an mp4 video of other codes such as the leapfrog or other numerical methods

    1. `Numerical_Intergrators.pu` :
        i  . SymplecticEulerStep(...)
        ii . LeapFrogStep(...)
        iii. YoshidaIntegratorStep(...)


Lint to create a better code