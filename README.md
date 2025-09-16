# TNNLS_Hybrid-Transfer-Active-Learning

This repository provides the codes required to generate the results for paper: "Hybrid Transfer Active Learning for Multi-Stream Processes with Within-Process and Cross-Process Correlation Modeling and Online Updating". 

To successfully run the code, you should download the nlopt package from the website "https://nlopt.readthedocs.io/en/latest/" and configure it for MATLAB usage following the installation manual for MATLAB: "https://nlopt.readthedocs.io/en/latest/NLopt_Installation/#matlab".

The codes beginning with "Figure" are utilized to generate the corresponding figures for different experiments, which will run the corresponding codes begin with "GenerateResults". For example, "FigTrig" is utilized to generate plots for the trigonometric signals, i.e., simulation one in the manuscript, within which the code "GenerateResultsTrigRandom" and "GenerateResultsTrigAL" are implemented to generate results for different methods under random and active learning scenarios. Also note that running each code beginning with "Figure" requires the running of NN results as well.  For the results of NN, please run the corresponding python codes.

The codes for 2D inputs are within the "Code2DInputs" folder.
