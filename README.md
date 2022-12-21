# Introduction
TRS is a efficient lattice-based traceable ring signature scheme. 

We give this code to test the computational overhead of TRS for the security parameter λ = 128.
## User Guide
Main_TRS is a main function, which calls the other five algorithms to realize the simulation experiment of the proposed scheme. Therefore, after downloading them, please run Main_TRS directly in MATLAB. Considering stability, for our core algorithm, we test 1000 times in the code. For more information, please see the comments in the various parts of the code.
### 1. Concrete steps
```
    Step 1: Run Main_TRS.
    Step 2: In the command line window, enter the ring size as prompted.
    Step 3: In the command line window, follow the prompts to enter the message to be signed.
```
### 2. Output
The output has three parts:

    1) The running time of key generation 1000 times.
    
    2) The running time of signature generation 1000 times.
    
    3）The judgment result of traceable ring signature validity and the running time of signature verification 1000 times.
