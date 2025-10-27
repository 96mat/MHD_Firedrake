An MHD solver based on [[1]](https://arxiv.org/pdf/2012.04122) and then stabilized using an augmented lagrangian term [[2]](https://arxiv.org/pdf/1706.02648) and small diffusive terms
needed to relax the stiffness of the monolithic problem (two-way coupled). Then an addition of the Hall's term

1. Code can be run in parallel    ```$ mpiexec -n 10 python First_test.py```
2. Advised mesh resolution ```128 x 128```
