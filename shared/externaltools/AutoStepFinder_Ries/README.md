This repository contains code associated with the paper:

*AutoStepfinder: a fast and automated step detection method for single-molecule analysis*

Luuk Loeff[1,2,3], Jacob W. J. Kerssemakers[1,3], Chirlmin Joo1 [4], Cees Dekker [1,4]

1 Kavli Institute of Nanoscience and Department of Bionanoscience, Delft University of Technology, Delft, The Netherlands
2 Present address: Department of Biochemistry, University of Zurich, Zurich, Switzerland
3 Equal contribution
4 Correspondence: C.Joo@tudelft.nl (CJ), C.Dekker@tudelft.nl (CD)

SUMMARY
Single-molecule techniques allow the visualization of the molecular dynamics of nucleic acids and proteins with high spatio-temporal resolution, for example a motor protein stepping along DNA. Valuable kinetic information of biomolecules can be obtained when the discrete states within single-molecule time trajectories are determined. Here we present a fast, automated, and bias-free step-detection method, AutoStepfinder, that we developed to determine steps in large datasets without requiring prior knowledge on the noise contributions, distribution, and location of steps. The analysis is based on a series of partition events that minimize the difference between the data and the fit. A dual-pass strategy determines the optimal fit and allows AutoStepfinder to detect steps of a wide variety of sizes. We demonstrate successful step detection for a broad variety of experimental traces. The user-friendly interface and the automated detection of AutoStepfinder provides a robust analysis procedure that enables anyone without programming knowledge to generate step fits and informative plots in less than an hour. 

HIGHLIGHTS
•	Fast, automated, and bias-free detection of steps within single-molecule trajectories.
•	Robust step detection without any prior knowledge on the data.
•	A dual-pass strategy for the detection of steps over a wide variety of scales.
•	A user-friendly interface for a simplified step fitting procedure.

REPOSITORY CONTENTS
•	code to analyze steps: 'AutoStepFinder'
•	code for cleaning data from corrupted points (such as Inf): 'DataDuster'
•	code for generating test traces : 'StepMaker'
•	zip file containing test traces

FOR USERS
•	an elaborate manual is included with Supplemental information in the paper
•	please cite the paper when using this code
•	for questions regarding this code, please contact:
	Dr. Jacob Kerssemakers (j.w.j.kerssemakers@tudelft.nl) or
	Dr. Luuk Loeff (l.loeff@bioc.uzh.ch)
