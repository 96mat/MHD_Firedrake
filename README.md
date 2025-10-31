# MHD Solver: Resistive and Hall preserving both divB and divU

An MHD solver based on [[1]](https://arxiv.org/pdf/2012.04122) and then stabilized using an augmented lagrangian term [[2]](https://arxiv.org/pdf/1706.02648) and small diffusive terms
needed to relax the stiffness of the monolithic problem (two-way coupled). Then an addition of the Hall's term

1. Code can be run in parallel    ```$ mpiexec -n 10 python First_test.py```
2. Advised mesh resolution ```128 x 128```

Choosing the right conforming DeRham complex, compatible spaces this solver was able to reproduce the classical cavity-Lid driven cavity problem. Some extracts (segmentations of the magnetic field is due to N1curl finite elements that are discontinuous by construction):

https://github.com/user-attachments/assets/348ec988-5479-4296-a5ac-901d04263663


https://github.com/user-attachments/assets/09173d8d-58f5-47a7-bdc3-b72a73e8f650




