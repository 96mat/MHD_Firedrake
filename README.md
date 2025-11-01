# Incompressible/Compressible MHD Solver: Resistive and Hall preserving both divB and divU

An MHD solver based on [[1]](https://arxiv.org/pdf/2012.04122) and then stabilised using an augmented lagrangian term [[2]](https://arxiv.org/pdf/1706.02648) and small diffusive terms
needed to relax the stiffness of the monolithic problem (two-way coupled). Then an addition of the Hall's term

1. Install ```Firedrake``` locally from [here](https://www.firedrakeproject.org/) or run it in ```GoogleColab``` attaining to the [following](https://github.com/firedrakeproject/firedrake/discussions/3302) procedure
2. Code can be run in parallel    ```$ mpiexec -n 10 python First_test.py```, through local installation
3. Advised mesh resolution ```128 x 128```
4. Choosing your finite element space combination depends on the boundary conditions you want to solve for
   
### TO DO:
1. Implement the compressibility i.e. a transport equation for the density that won't be constant anymore
2. Possible coupling with other physics, like the Temperature and diffusion-reaction-transport of species
3. To precondition the fluid and magnetic blocks separately, enhancing the possible parameters' range for which the solver would be suited to "solve-for"
   
# Some Incompressible extracts
<img width="883" height="369" alt="image" src="https://github.com/user-attachments/assets/930d55f6-8b83-496a-8ef2-6fdf0bdcf179" />


Choosing the right conforming DeRham complex, compatible spaces, this solver was able to reproduce the classical lid-driven  cavity problem. Some extracts (segmentations of the magnetic field are due to N1curl, Raviart-Thomas finite elements that are discontinuous by construction):

https://github.com/user-attachments/assets/348ec988-5479-4296-a5ac-901d04263663


https://github.com/user-attachments/assets/ef9d8c3a-c304-4be0-87a2-dbae2f74b8da


https://github.com/user-attachments/assets/09173d8d-58f5-47a7-bdc3-b72a73e8f650


https://github.com/user-attachments/assets/3ae82276-631a-43b7-bf8a-7504473aa9e8


