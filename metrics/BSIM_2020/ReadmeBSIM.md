# BSIM2020 Model
The BSIM2020 Model uses a blind binaural processing stage for simulating 
human binaural processing.
It requires a two-channel (left/right) binaural signal, which contains 
the mixture of speech and noise. The model's output is a structure 
containing the binaural enhanced signal. The model applies two processing 
strategies of the binaural system, which run in parallel. 
a) level minimization, which attenuates the dominant signal (i.e. noise at 
negative SNRs)
b) level maximization, which amplifies the dominant signal (i.e. the target
at positive SNRs).
The best of both strategies is selected blindly by employing a modified version
of the SRMR model, which is applied on a channel-by-channel basis.  
The same SRMR model is applied to select the better ear.
Please note that the modified SRMR model 
either considers the whole signal as a single time frame (long
term model) or uses the same time constants as the EC processing and BE
processing if the short time model is used.
The time constant for short time EC processing is 300ms [with 50% overlap,
leading to an effective window of 150ms (see Hauth and Brand (2018) for further
details] and the time constant for short time Better-Ear processing is 20 ms. 

The SRMR measure is described here:
(Tiago H. Falk, Chenxi Zheng, and Way-Yip Chan. A Non-Intrusive Quality and 
Intelligibility Measure of Reverberant and Dereverberated Speech,
IEEE Trans Audio Speech Lang Process, Vol. 18, No. 7, pp. 1766-1774, 
Sept. 2010. 
[doi:10.1109/TASL.2010.2052247](http://dx.doi.org/10.1109/TASL.2010.2052247))

The SRMR measure can be downloaded here:
https://github.com/MuSAELab/SRMRToolbox

This model also uses the Large Time-Frequency Analysis Toolbox:
https://ltfat.github.io/

# Example:
An Example of how to use the model can be found in "DemoMain.m".

A set of Olsa sentences is provided (in the folder "Anechoic"). The recordings
were kindly provided by Sabine Hochmuth. 
Please note that they are not taken from the original Olsa corpus. Therefore,
predictions might slightly differ from those with the original Olsa sentences.
 
If you want to use original Olsa sentences, please
contact HörTech gGmbH (https://www.hoertech.de/en/home-ht.html).

The BSIM2020 can be used in two different modes:
1) Processing of mixed signals only:
out_struct = BSIM20('RequiredSignal',mixed_input,'model_params',model_param);
Here, the mixed input signals are blindly processed by the binaural processing
stage. The outputis a binaurally enhanced signal.

Parameters:
'Required Signals'
mixed_input is a two-column matrix, where the left coulumn is associated to
the left ear and the right columns is associated to the right ear. 

In the output structure out_struct, the following signals are available:
1) Blindly reconstructed binaurally processed signal
2) Output of the minimization strategy + better ear
3) Output of the maximization strategy + better ear
4) Left ear only -> same signal as input
5) Right ear only -> same signal as input

2) An additional set of signal can be provided to the model.
They are denoted as 
'Optional Signals'
Optional signals, which are processed in the same way as the mixed signals
(sometimes this processing is referred to as "shadow filtering").
In this case, the BSIM20 is invoked with:
out_struct = BSIM20('RequiredSignal',mixed_input,'OptionalSignal',OptionalSignals,'model_params',model_param));

Optional Signals can contain an arbitrary number of signals. The only requirement
is that they are put pair-wise into the matrix. 
[L R L R ---L R
 E I E I --- .
 F G F G --- .
 T H T H --- .
 1 T 2 T ---N
   1   2     N ]

# Speech Intelligibility Back-end
This software package comes with a modified version of 
the Speech Intelligibility Index (SII)in order to demonstrate, how a 
non-blind back end can be used with this model.
The SII requires frequency-band specific levels of target speech and 
interfering signals in order to work correctly. In this release, it is 
assumed that the first optional signal is speech, the second optional signal is noise. 
If optional signals are present, the output structure also contains the frequency
specific levels. The user is responsible that he/she knows, which optional signal
is the target signal.
The original matlab SII implementation can be found here: http://sii.to/programs.html


# References

If you use this model in your research, please cite the reference below:

Hauth, C.F., Berning, S. C., Kollmeier, B., Brand, T., 2020. 
Modelling binaural unmasking of speech using a blind binaural processing stage
Trends Hear.
# Contact:
If you have any questions feel free to contact either me 
christopher.hauth@uni-oldenburg.de
or my colleague Dr. Thomas Brand
thomas.brand@uni-oldenburg.de

%--------------------Licence ---------------------------------------------
% Copyright (c) <2020> Christopher F. Hauth
% Dept. Medical Physics and Acoustics
% Carl von Ossietzky University Oldenburg 
% Permission is hereby granted, free of charge, to any person obtaining 
% a copy of this software and associated documentation files 
% (the "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to 
% permit persons to whom the Software is furnished to do so, subject 
% to the following conditions:
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% END OF FILE