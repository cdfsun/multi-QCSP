# multi-QCSP
Benders decomposition for quay crane scheduling problem

File "multi_QCSP_model.cpp": the C++ code of QCSP model, implemented through CPLEX.

File "multi_QCSP_Benders.cpp": the C++ code of Benders approach, including the lower bound provided by the unidirectional formulation, implemented through CPLEX.

File "test": the file of test data from real practice. All of the above instances are named by "n-b-q", where n is the number of tasks, b is the number of bays, and q is the mumber of QCs. In the txt file corresponding to each instance, the first "[]" records the number of tasks, the number of bays, the number of precedence pairs, the number of QCs, the travel time, safety margins. The next two "[]"s record the processing times and bay locations of each task, respectively. The following two "[]"s record the ready times, initial bay locations of each QC, respectively. The rest "[]"s record the detailed task precedence pairs.
