%%% readme.txt for AN model %%%

This is the BEZ2018a version of the code for auditory periphery model from
the Carney, Bruce and Zilany labs.

This release implements the version of the model described in:

	Bruce, I.C., Erfani, Y., and Zilany, M.S.A. (2018). "A Phenomenological
	model of the synapse between the inner hair cell and auditory nerve: 
	Implications of limited neurotransmitter release sites," Hearing Research 360:40-54.
	(Special Issue on "Computational Models in Hearing".)

with the synapse modifications described in:

	Bruce, I., Buller, A., and Zilany, M. "Modeling of auditory nerve fiber input/output functions
	near threshold," Acoustics 2023, Sydney, Australia, December 2023.

Please cite these two publications if you publish any research results obtained with
this code or any modified versions of this code.

*** Change History ***

See the file changelog.txt

*** Instructions ***

The Matlab and C code included with this distribution is designed to be
compiled as a Matlab MEX file, i.e., the compiled model MEX function will run
as if it were a Matlab function.  The code can be compiled within Matlab using
the function:

    mexANmodel.m

Note that it is also possible to compile and run the code in the open-source
(i.e, free!) software Octave just as it is done in Matlab, although it will
typically run more slowly in Octave.  It may be necessary to install the
"signal" and "control" packages within Octave before running the AN model code.

The code is implemented in two parts. The first function, model_IHC_BEZ2018a,
takes the acousic signal as input and gives the IHC's relative transmembrane
potential (Vihc) as the output.  The second function, model_Synapse_BEZ2018a,
takes Vihc as input and outputs the PSTH (or a spike train for a single stimulus
presentation).  There are also a number of additional optional outputs from
this second function - see the help information for further details. For
instructions on how to run the code, the following commands can be run at
the Matlab command prompt:

    help model_IHC_BEZ2018a

    help model_Synapse_BEZ2018a

We have also included:-

1. a Matlab function "generateANpopulation.m" for generating a population
   of AN fibers with the statistics for spont rate and absolute & relative
   refractory periods as described in the journal article,

2. a function "fitaudiogram2.m" for estimating the parameters for outer and
   inner hair cell impairment, Cohc and Cihc, respectively, for a given
   audiogram, and

3. a number of scripts with filenames beginning with "test..." that run
   simulations and generate figures similar to a subset of those from the
   journal article and supplementary material.

ACKNOWLEDGMENTS

Earlier contributions to the model code made by Xuedong Zhang,
Michael G. Heinz, Qin Tan, Paul C. Nelson, and Laurel H. Carney.

%%% Ian C. Bruce (ibruce@ieee.org), Muhammad S. A. Zilany (msazilany@gmail.com)
           - December 2023 %%%
