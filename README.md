# Incompressible/Compressible MHD Solver: Resistive and Hall preserving both divB and divU

A 2D MHD solver based on [[1]](https://arxiv.org/pdf/2012.04122) and then stabilised using an augmented lagrangian term [[2]](https://arxiv.org/pdf/1706.02648) and small diffusive terms
needed to relax the stiffness of the monolithic problem (two-way coupled). Then an addition of the Hall's term

1. Install ```Firedrake``` locally from [here](https://www.firedrakeproject.org/) or run it in ```GoogleColab``` attaining to the [following](https://github.com/firedrakeproject/firedrake/discussions/3302) procedure
2. Code can be run in parallel    ```$ mpiexec -n 10 python First_test.py```, through local installation
3. Advised mesh resolution ```128 x 128```
4. Choosing your finite element space combination depends on the boundary conditions you want to solve for
   
### TO DO:
1. Implement the compressibility i.e. a transport equation for the density that won't be constant anymore :negative_squared_cross_mark:
2. Possible coupling with other physics, like the Temperature and diffusion-reaction-transport of species :negative_squared_cross_mark:
3. To precondition the fluid and magnetic blocks separately, enhancing the possible parameters' range for which the solver would be suited to "solve-for" and errors :white_check_mark:
4. Solve the problem in its Fully-Non-Linear statement :white_check_mark:
   
# Some extracts from the incompressible formulation 
<img width="883" height="369" alt="image" src="https://github.com/user-attachments/assets/930d55f6-8b83-496a-8ef2-6fdf0bdcf179" />


Choosing the right conforming DeRham complex, compatible spaces, together with the augmented Lagrangian term, this solver was able to preserve $\nabla\cdot\boldsymbol{B}$ and reproduce the classical lid-driven cavity problem. Some extracts (segmentations of the magnetic field are due to N1curl, Raviart-Thomas finite elements that are discontinuous by construction):
## :one: $H^1(\Omega)/H^1(curl;\Omega)$ and $H^1(div;\Omega)/H^1(div;\Omega)$ and MHD: only Resistive
These two combinations of finite elements are not DeRham compatible; therefore, a total conservation of $\nabla\cdot B$ and $\nabla \cdot u$ is not expected and unfortunately verified (in factions of cents), but the need to use these spaces comes from their extreme flexibility when it comes to imposing boundary conditions (strongly), especially when using the classical $H^1(\Omega)$ continous finite elements

https://github.com/user-attachments/assets/348ec988-5479-4296-a5ac-901d04263663


https://github.com/user-attachments/assets/ef9d8c3a-c304-4be0-87a2-dbae2f74b8da


https://github.com/user-attachments/assets/09173d8d-58f5-47a7-bdc3-b72a73e8f650

## :two: $H^1(div;\Omega)/H^1(curl;\Omega)$ compatible spaces and MHD: comprising the Hall's Term
Here, employing the fully DeRham complex, i.e. H-div and H-curl combo, allows us to preserve $\nabla\cdot B$ and $\nabla\cdot u$ up to machine precision, but, for instance, the obiquitous no-slip BC condition doesn't exist by default in Raviart-Thomas spaces and has to be implemented following the Nitsche weakly imposition typical of Discontinuous Galerkin elements [[3]](https://math.okstate.edu/people/yqwang/publications/divfree-hdiv-stokes-v3.pdf) -that only works for small $Re$ i.e. high $\nu$ values

https://github.com/user-attachments/assets/3ae82276-631a-43b7-bf8a-7504473aa9e8


https://github.com/user-attachments/assets/d2acd7eb-20a6-4e7d-8b5b-a8bdf42eae42

The div-free cavity lid problem in RT/Nedèlec finite elements, at $Re=1$ 



https://github.com/user-attachments/assets/77c0fffc-a954-4662-b662-58d4f872b223

