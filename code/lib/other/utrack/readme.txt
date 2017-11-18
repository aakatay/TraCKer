
run in the parent folder

-crosb004 : (combined run)
calls Scott's code LongMultiMovieSplitAnalysis
change: 
exp_name = '150423cell15';
sectionsize = 200;
framegap = [2]; 


-LongMultiMovieSplitAnalysis: 
calls CME in a for loop for multi data analysis
change: 
cmeAnalysis('Parameters', [1.45, 100, 16] ...
[NA magnification pixel size]

-prob004 : (post run _ before 004)
combines the output of the crosb004 from the Section folders
