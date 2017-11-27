The codes in this zip file allow you to estimate a model in WTP space with a flexible distribution of WTPs and the price/scale coefficient. You need to have Matlab and its optimization toolbox on your machine in order to run the codes. Runs will be much faster if you also have Matlab's parallel computing toolbox on your machine along with a fast gpu; but the codes will run without them.

The method is explained in this paper:
K. Train, 2015, "Mixed Logit with a Flexible Mixing Distribution," working paper, available here:
http://elsa.berkeley.edu/~train/flexible.pdf

Models in WTP space are discussed in this paper:
K. Train and M. Weeks, 2005, "Discrete Choice Models in Preference Space and Willingness-to-Pay Space," Ch. 1, pp. 1-17, in Applications of Simulation Methods in Environmental Resource Economics, A. Alberini and R. Scarpa, eds., Springer Publisher: Dordrecht, The Netherlands, available here
http://eml.berkeley.edu/~train/trainweeks.pdf


The zip file contains numerous matlab .m files. The file FlexibleWtp.m is the only one you need to change in order to run your models. The other .m files are called directly or indirectly by FlexibleWtp.m. You can put all the files into one folder, and run FlexibleWtp.m from that folder. Or you can put all of the .m files except FlexibleWtp.m into a different folder and add a command to FlexibleWtp.m that adds that folder as a path.  

To run your model, specify the terms in code FlexibleWtp.m. The code gives instructions for each term that you need to specify. After you have changed FlexibleWtp.m for your model, execute the file. It will perform the estimation for the model and data that you specify.

The zip file contains a sample dataset called videodata100.mat, which contains the choices of the first 100 respondents in the survey data used in the paper. FlexibleWtp.m is currently specified to run a model on this dataset. The output for this run is shown in MyResults_KT.out. Before trying to run your own models, run FlexibleWtp.m as is and check to see that the output you get matches MyResults_KT.out. (Note: the code is currently set with YesGPU=0 which tells Matlab not to use the parallel processing toolbox, because some users might not have this toolbox. If you have Matlab's parallel processing toolbox, then, after running the example model and checking the output, be sure to change YesGPU to 1 to get faster run times for your own models. Also, NReps is set to only 4 for bootstrapping, so that the example model will run quickly for you. Bootstrapping usually uses more than 4 resamples.)

The code does not currently allow estimation of models in preference space. I might provide that option at some later date.

If you find any bugs in the program, please email me about them.  

If you want to be notified of any changes in the code, including fixes to bugs, then please email me. I'll keep a list of users and notify everyone on the list when/if new versions of the software are available.

Best wishes,
Kenneth Train
train@econ.berkeley.edu  
