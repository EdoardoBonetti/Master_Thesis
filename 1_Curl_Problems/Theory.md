# Decomposition of a vector Field
In calculus it is shown how any smooth field in $\mathbb{R}^3$ can be divided into two parts:
$$
\begin{align}
C^{\infty}(\Omega; \mathbb{R}^3) = \nabla C_0^\infty(\Omega) \otimes \nabla\times C^\infty(\Omega;\mathbb{R}^3 ) 
\end{align}
$$
to be more precise for $u \in C^{\infty}(\Omega; \mathbb{R}^3) $ there will be decomposition $(\phi, \psi) \in C_0^{\infty}(\Omega) \times C^{\infty}(\Omega; \mathbb{R}^3) $ called **potentials** such that :
$$
\begin{cases}
u = \nabla \phi + \nabla \times \psi \\
\nabla \cdot \psi = 0 \\
\end{cases}
$$
In particular the construction follows this idea:
1. We solve the poisson problem :
$$
\begin{cases}
\Delta \phi = - \nabla\cdot u  \\
\partial_n \phi = -u
\end{cases}
$$  
From the theory we know that There is a solution to this problem.

2. By construction know that now the function $u+\nabla \phi$ is a divergence free function since $\nabla \cdot (u+\nabla \phi) =\nabla \cdot u - \nabla \cdot u =0 $ and on the boundary vanishes (for the same reason). Thanks to the Helmholtz decomposition theorem the kernel of the divergence is equal to the image of the curl operator, therefore there will be a three dimensional vector field $\psi$ such that $\nabla \times \psi = u+\nabla \phi $. 

The reader may notice that the decomposition is not unique, infact given $(\phi, \psi)$ and $(\phi', \psi')$ potentials we have that the first 2 differ of a laplace field, i.e. the $\nabla f$ where $\Delta f =0$ with vanishing Neumann b.c. . Furthermore, by construction, $\psi$ and $\psi'$ differ too.

There exists a method to obtain a unique decomposition setting certain boundary conditions such as:
1. Homogeneous dirichlet b.c. on the boundary of the domain for the scalar potential.
2. Homogeous neumann b.c. on the boundary of the domain for the vector potential.

In this way one can obtain a unique decomposition thanks to the inf-sup theorem, i.e. the inf of the solutio of a Poisson problem  is less or equal to the inf of the boundary values, simialrly fo the sup part. 
In particular sice they are zero the difference between 2 possible solution must be zero $\implies$ uniqueness.

A similar approach can be used to determine the more general decomposition of a vector field in $\mathbb{L}^2$ , in particular the approach is presented in the notes on Numerical Methods for Maxwell's Equations, by Joachim Schoeberl, subsection 2.3.2 pp. 24-27 : https://www.asc.tuwien.ac.at/~schoeberl/wiki/lva/notes/numpde.pdf .
The reader may notice that there will be different possibilities to decompose a vector field in $\mathbb{L}^2$ . Therefore one may pick the decomposition that is most appropriate for the problem.

Here I will state the decomposition as it is reported in the notes and I will pick one possible decomposition to implement.

Let $q\in [ \mathbb{L}^2(\Omega) ] ^3$ There exists a decomposition
$$
\begin{align}
q = \nabla \phi + \nabla \times \psi 
\end{align}
$$
The following choices for the function $(\phi,\psi)$ are possible (they are mutually exclusive), and are bounded by $\| q \|_{\mathbb{L}^2}$ in their respective norm:
1. $\phi \in H^1 \text{ and } \psi \in [H^1]^3  \text{ s.t. div}\psi = 0 $

2. $\phi \in H^1 \text{ and } \psi \in [H_0^1]^3$

3. $\phi \in H^1 \text{ and } \psi \in H_0(\text{curl})  \text{ s.t. div}\psi = 0$

4. $\phi \in H_0^1 \text{ and } \psi \in [H^1]^3 \text{ s.t. div }\psi = 0$

5. $\phi \in H_0^1 \text{ and } \psi \in H(\text{curl }) \text{ s.t. div }\psi = 0 \text{ and } tr_n \psi =0$

In case one is wondering on how to implement the above I strongly recommend to read the notes mentioned above. The proofs are all constructive with the only exeption where the fourier thansform is used, in that case we use a more numerical approach that consists of a mixed formulation to project the problem into the corrispective space.

In particular the one implementes in the jupyter notebook is number (3), but it does not follow the constructive proof if the notes , it is a simple projection onto subspaces of $[L_2]^3$.

