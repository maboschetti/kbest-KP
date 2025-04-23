++++++++++++  Code and executable  ++++++++++++ 

Folder "Code" includes the files:

kBest-KP-Library.cpp
kBest-KP-Library.h
Main-KP-Test.cpp
Utilities.cpp
Utilities.h

By building the code one should obtain an executable. Call it "KP-kbestLib-v3.exe" 
and put it in the folder "Tests" and "TestsLargeInstances".

++++++++++++  Algorithm Name  ++++++++++++ 

We use abbreviations for the algorithm types:

ALGORITHM KPFW_kbest_L is called A4
ALGORITHM KPBW_kbest_L is called A7 
ALGORITHM KPBW_kbest_LDCF(M) is called A10

++++++++++++  Instances  ++++++++++++

Folder "Instances":

Includes the instance files (both normal size and large size)
Each file name reports: n_w_M_D_run.txt

++++++++++++  Tests on normal instances ++++++++++++  

Folder "Tests":

Includes all the files .bat for the tests.

Must include the executable: KP-kbestLib-v3.exe

To run the tests one should run file tests.bat 

Includes a folder "Results", where the executable will print results on .cvs files.
A file for algorithm type, k value, and number of the run (three runs).

The output file for the algorithm A4 shows these values:
-the name of the instance
-the algorithm type
-the id of the instance
-the id of the instance randomly generated
-the number of the run (0, 1, 2)
-n
-W
-k
-M
-D
-the best solution
-the kth solution (or the last solution found if k solutions do not exist)
-the generation computing time
-the total computing time

The output files for the algorithms A7 and A10 show these values:
-the name of the instance
-the algorithm type
-the id of the instance
-the id of the instance randomly generated
-the number of the run (0, 1, 2)
-n
-W
-k
-M
-D
-the best solution
-the kth solution (or the last solution found if k solutions do not exist)
-the average k
-the number of states used
-the % of states used
-the number of states generated
-the % of states used
-the average dimension of k
-the % dimension of k
-the generation computing time
-the total computing time

++++++++++++  Tests on large instances ++++++++++++  

Folder "TestsLargeInstances":

Includes all the files .bat for the tests.

Must include the executable: KP-kbestLib-v3.exe

To run the tests one should run file "testsA10_large.bat" 

Includes a folder "Results", where the executable will print results on .cvs files.
A file for algorithm type, k value, and number of the run (three runs).

The output files for the algorithm A10 shows these values:
-the name of the instance
-the algorithm type
-the id of the instance
-the id of the instance randomly generated
-the number of the run (0, 1, 2)
-n
-W
-k
-M
-D
-the best solution
-the kth solution (or the last solution found if k solutions do not exist)
-the average k
-the number of states used
-the % of states used
-the number of states generated
-the % of states used
-the average dimension of k
-the % dimension of k
-the generation computing time
-the total computing time

