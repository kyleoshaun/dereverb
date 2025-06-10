% model_Synapse_BEZ2018a - Bruce, Erfani & Zilany (2018) Auditory Nerve Model (Dec 2023 update)
%
%     psth = model_Synapse_BEZ2018a(vihc,CF,nrep,dt,noiseType,implnt,spont,tabs,trel,expliketype);
%
% psth is the peri-stimulus time histogram (PSTH) (or a spike train if nrep = 1)
%
% vihc is the inner hair cell (IHC) relative transmembrane potential (in volts)
% CF is the characteristic frequency of the fiber in Hz
% nrep is the number of repetitions for the psth
% dt is the binsize in seconds, i.e., the reciprocal of the sampling rate (see instructions below)
% noiseType is for "variable" or "fixed (frozen)" fGn: 1 for variable fGn and 0 for fixed (frozen) fGn
% implnt is for "approxiate" or "actual" implementation of the power-law functions: "0" for approx. and "1" for actual implementation
% spont is the spontaneous firing rate in /s
% tabs is the absolute refractory period in s
% trel is the baselines mean relative refractory period in s
% expliketype sets the type of exp-like function in IHC -> synapse mapping: 1 for shifted softplus (preferred); 0 for no expontential-like function; 2 for shifted exponential; 3 for shifted Boltmann
% 
% For example,
%
%    psth = model_Synapse_BEZ2018a(vihc,1e3,10,1/100e3,1,0,50,0.7,0.6,1); **requires 10 input arguments
%
% models a fiber with a CF of 1 kHz, a spontaneous rate of 50 /s, absolute and baseline relative
% refractory periods of 0.7 and 0.6 s, for 10 repetitions and a sampling rate of 100 kHz, and
% with variable fractional Gaussian noise, the approximate implementation of the power-law functions
% in the synapse model, and the shifted softplus function in the IHC -> synapse mapping
%
% OPTIONAL OUTPUT VARIABLES:
%
%     [psth,meanrate,varrate,synout,trd_vector,trel_vector] = model_Synapse_BEZ2018a(vihc,CF,nrep,dt,noiseType,implnt,spont,tabs,trel,expliketype);
%
% meanrate is the analytical estimate of the mean firing rate in /s for each time bin
% varrate is the analytical estimate of the variance in firing rate in /s for each time bin
% synout is the synapse output rate in /s  for each time bin (before the effects of redocking are considered)
% trd_vector is vector of the mean redocking time in s for each time bin
% trel_vector is a vector of the mean relative refractor period in s for each time bin
%
% NOTE ON SAMPLING RATE:-
% Since version 2 of the code, it is possible to run the model at a range
% of sampling rates between 100 kHz and 500 kHz.
% It is recommended to run the model at 100 kHz for CFs up to 20 kHz, and
% at 200 kHz for CFs> 20 kHz to 40 kHz.