function out_struct = BSIM20(varargin)
% Usage: out_stuct = BSIM20('RequiredSignal',mixed_input,'OptionalSignal',...
% [OptionalSignal1 OptionalSignal2]...,OptionalSignalN,'model_params',model_params)
%
% Input:
% 'RequiredSignal' is the mixture of speech and noise, which is a two
% channel (left,right) matrix (mixed_input)
% 'OptionalSignal' - One or more optional two channel signals, which are
% processed in the same way as the 'RequiredSignal'. They have to be
% arranged in a single matrix such that always two columns correspond to
% the left and right ear signals [left1 right1 left2 right2 ... leftN rightN]
% (e.g. if you need clean speech and noise for your back-end)
% 'model_params' - Structure containing model parameters like frequency
% range etc.
%
% Output: 
% out_struct -  Structure containing the processed signals, in the order [Sel Min Max]
% e.g. out_struct.signals.Mixsigsynfin contains the EC/BE processed mixed
% signal, where the best strategy was chosen blindly
% The SII used in the study by Hauth et al. (2020) 
% requires band specific levels. If optional signals are present,
% the levels of all optional signals after EC processing and better ear
% selection are part of the output structure.

% Authors: Christopher F. Hauth <christopher.hauth@uni-oldenburg.de>
%          Dr. Thomas Brand     <thomas.brand@uni-oldenburg.de> 
% Date: 22.10.2020
%--------------------------------------------------------------------------
%% Input check the variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_flags  = zeros(4,1);
for kk=1:2:length(varargin)-1
    switch varargin{kk}
        case 'RequiredSignal'
            % Contains the mixed signals:
            MixedSig = varargin{kk+1};
            param_flags(1) = 1;
        case 'OptionalSignal'
            % Contains the optional signals:
            OptSigs = varargin{kk+1};
            iNumofOptsigs = size(OptSigs,2)/2;
            param_flags(2) = 1;
        case 'model_params'
            % contains model parameter:
            model_param = varargin{kk+1};
            fs = model_param.fs;
            param_flags(3) = 1;
            frange = [model_param.fmin  model_param.fmax];
            BMtype = model_param.Filterbank;
            batch_flag = model_param.long_term;
            bin_inaccuracy  = model_param.bin_err;
            ERB_factor = model_param.ERB_factor;
        otherwise
            return
    end
end
% if one of the inputs is missing, replace either by default value or throw
% a warning or error message with a description what to do:
flag_idcs = find(param_flags==0, 1);
if isempty(flag_idcs)
    disp('Parameter check complete! Processing starts....')
else 
    if param_flags(1)==0
        error('Mixed signals are required as input!')
    elseif param_flags(2)==0
        %disp('Only mixed signals are processed, no optional signals are provided.'); Removed by Kyle O
    elseif param_flags(3)==0
        frange = [125 8000];
        disp('Missing Parameter: Frequency range is 150 Hz to 8000 Hz.');
        BMtype = 'GT';
        disp('Missing Parameter: Gammatone filterbank is used.');
        batch_flag = 1;
        disp('Missing Parameter: Batch processing is used.');
        bin_inaccuracy = 1;
        disp('Missing Parameter: Binaural processing inaccuracies are used.');
        ERB_factor = 1; 
        disp('Missing Parameter: ERB factor of 1 is used.');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_model            = 44100;% sampling rate of model
if frange(1)==frange(2)
    ftarget         = frange(1);
else 
    ftarget         = 500; %[Hz]    
end
 % resmaple input if necessary
if fs_model~=fs
    MixedSig        = resample(MixedSig,fs_model,fs);
    if exist('OptSigs','Var')
        OptSigs     = resample(OptSigs,fs_model,fs);
    end
end
% define temporal windows for EC processing and better ear selection
% batch processing means that the whole signal is considered as a single time frame
if batch_flag 
    iBlockLen_EC    = length(MixedSig); 
    iBlockLen_BE    = length(MixedSig);
    iNumofBlocks_EC = 1;
    iNumofBlocks_BE = 1;
else
    % If short time processing is performed, EC process and better ear
    % processing operate on different time constants.
    % The EC process is conducted in
    % time frames of 300ms (according to Hauth&Brand, 2018)
    % Using a frame shift of 50% leads to an effective temporal window of
    % 150 ms.
    iBlockLen_EC    = round(0.300.*fs_model); 
    % Better ear processing is performed on time frames of 20ms
    % Using a frame shift of 50% leads to an effective temporal window of
    % 10 ms.
    iBlockLen_BE    = round(0.020.*fs_model);
    
    % make sure the block length is even
    if mod(iBlockLen_EC,2)
        iBlockLen_EC = iBlockLen_EC+1;
    end
    iNumofBlocks_EC = 2.*floor(length(MixedSig)./iBlockLen_EC)-1;
    iNumofBlocks_BE = 2.*floor(length(MixedSig)./iBlockLen_BE)-1;
    iBlockshift_EC  = floor(iBlockLen_EC/2); % 50 percent overlap
    iBlockshift_BE  = floor(iBlockLen_BE/2);
    window_EC       = hann(iBlockLen_EC,'periodic');
    window_BE       = hann(iBlockLen_BE,'periodic');
end
% Binaural Processing Inaccuracies defined by vom Hövel (1984)
% For details see supplementary material of Hauth et al. (2020) (not yet published)
if bin_inaccuracy
    % error in level 
    sigma_epsilon0 = 1.5; % [dB]
    alpha0         = 15;  % [dB]
    p              = 1.6; % no unit as it used as exponent
    
    % error in time
    sigmadelta0    = 65.*10^-6; % [sec]
    Delta0         = 1.6.*10^-3;% [sec] 
else
    sigmadelta0    = 0;% [sec]
    Delta0         = 1;
end
% Parameters for the gammatone filterbank
order              = 4;
filter_per_ERB     = 1;

%parameter defining the cut-off frequency of EC processing. Above this
%frequency, better-ear processing is assumed to play the dominant role.
f_cut = 1500; % 1500 Hz

% Gammatone Filterbank from described in Hohmann (2002)
if strcmp(BMtype,'GT')
    analyzer1 = gfb_analyzer_new(fs_model,frange(1),ftarget,frange(2),filter_per_ERB,order,ERB_factor);
     %% Apply the Gammatone filterbank
     % LEFT EAR
     [MixsigLbatch, ~] = gfb_analyzer_process(analyzer1, MixedSig(:,1)');
     fc = analyzer1.center_frequencies_hz;
     MixsigLbatch = real(MixsigLbatch)';
     SomeImagSig = imag(MixsigLbatch)';
     % RIGHT EAR
     [MixsigRbatch, ~] = gfb_analyzer_process(analyzer1, MixedSig(:,2)');
     MixsigRbatch = real(MixsigRbatch)';
     % If optional signals are present, process in the same way as mixed
     % signals
     if exist('OptSigs','var')
          for mm=1:iNumofOptsigs
                eval(sprintf('[OptSig%dLbatch, ~] = gfb_analyzer_process(analyzer1, OptSigs(:,2*(mm-1)+1));',mm));
                eval(sprintf('OptSig%dLbatch = transpose(real(OptSig%dLbatch));',mm,mm));
                eval(sprintf('[OptSig%dRbatch, ~] = gfb_analyzer_process(analyzer1, OptSigs(:,2*mm));',mm));
                eval(sprintf('OptSig%dRbatch = transpose(real(OptSig%dRbatch));',mm,mm));
          end
     end
end
% Allocate buffer for the final resynthesis of signals
% the ending fin denotes final representation of the signal
MixsigLfin    = zeros(size(MixsigLbatch));
MixsigRfin    = MixsigLfin;
Mixsigminfin  = MixsigLfin;
Mixsigmaxfin  = MixsigLfin;
Mixsigsynfin  = MixsigLfin;
% Allocate buffer for optional signals if present:
if exist('OptSigs','var')
    for ll=1:iNumofOptsigs
        eval(sprintf('OptSig%dLfin  = MixsigLfin;',ll));
        eval(sprintf('OptSig%dRfin  = MixsigLfin;',ll));
        eval(sprintf('OptSig%dminfin = MixsigLfin;',ll));
        eval(sprintf('OptSig%dmaxfin = MixsigLfin;',ll));
        eval(sprintf('OptSig%dsynfin = MixsigLfin;',ll));
    end
end
%% Binaural stage processing starts here:
% iterate through frequencies
for kk = 1:length(analyzer1.center_frequencies_hz)
    % Check, if frequency is in the range of frequencies where binaural
    % processing is assumed
    if analyzer1.center_frequencies_hz(kk)<=f_cut
        % Iterate through time blocks for EC processing
        for mm = 1:iNumofBlocks_EC  
            if iNumofBlocks_EC==1 % i.e. considering the signal as a single time frame
                % save signal portion 
                MixsigL = MixsigLbatch(:,kk);
                MixsigR = MixsigRbatch(:,kk);
                if exist('OptSigs','Var')
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dL = OptSig%dLbatch(:,kk);',ll,ll));
                        eval(sprintf('OptSig%dR = OptSig%dRbatch(:,kk);',ll,ll));
                    end
                end       
            else
                % Do block-wise processing 
                MixsigL = MixsigLbatch((mm-1)*iBlockshift_EC+1:...
                        iBlockLen_EC+(mm-1)*iBlockshift_EC,kk).*sqrt(window_EC);
                MixsigR = MixsigRbatch((mm-1)*iBlockshift_EC+1:...
                        iBlockLen_EC+(mm-1)*iBlockshift_EC,kk).*sqrt(window_EC);
                    % Do the same for optional signals if present
                    if exist('OptSigs','Var')
                        for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dL = OptSig%dLbatch((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk).*sqrt(window_EC);',ll,ll));
                        eval(sprintf('OptSig%dR = OptSig%dRbatch((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk).*sqrt(window_EC);',ll,ll));
                        end
                    end  
            end
            % Get levels of left and right ear signal
            rmsL = sqrt(mean(MixsigL.^2));
            rmsR = sqrt(mean(MixsigR.^2));
            % calculate mean between ears. This will be the target level
            % for the re-synthesis
            rms_target = mean([rmsL' rmsR'],2);
            
            % Level difference between left and right ear.
            alpha = 20.*log10(rmsL./rmsR);
            
            % Calculate level error (depending on the level difference)
            if bin_inaccuracy
                 sigma_epsilon = sigma_epsilon0.*(1+(abs(alpha)/alpha0).^p);
                 levErrorL = mean(10.^(sigma_epsilon.*randn(1,1)./20));
                 levErrorR = mean(10.^(sigma_epsilon.*randn(1,1)./20));
            else
                % if processing is assumed to be perfect:
                levErrorL = 1;
                levErrorR = 1;
            end
            % Calculate factors for level equalization:
            facR = sqrt((rmsL./rmsR).*levErrorR*levErrorL);
            facL = 1./facR;
           
            %% Apply equalization in Level (First step of EC Process)
            MixsigLproc = facL.*MixsigL;
            MixsigRproc = facR.*MixsigR;
            % Apply the same equalization to the Optional Signals:
            if exist('OptSigs','Var')
                for ll=1:iNumofOptsigs 
                    eval(sprintf('OptSig%dLproc=facL.*OptSig%dL;',ll,ll));
                    eval(sprintf('OptSig%dRproc=facR.*OptSig%dR;',ll,ll));
                end
            end
    %% Calculate delays to maximize and minimize the output of the EC
     % mechanism
           % run EC processing for the mixed signals:
           [Mixsigmin, Mixsigmax, ECparam4OptSigs] = fftcon(MixsigRproc,MixsigLproc,...
           analyzer1.center_frequencies_hz(kk),fs_model,sigmadelta0,Delta0,bin_inaccuracy);     
           
           if exist('OptSigs','Var')
               for ll=1:iNumofOptsigs
                    eval(sprintf('[OptSig%dmin,OptSig%dmax] = ECprocess4OptSigs(OptSig%dRproc,OptSig%dLproc,ECparam4OptSigs,fs_model,bin_inaccuracy);',ll,ll,ll,ll));
               end
           end
 
        rms_max   = rms(Mixsigmax);
        Mixsigmax = rms_target./rms_max.*Mixsigmax;
        
        rms_min   = rms(Mixsigmin);
        Mixsigmin = rms_target./rms_min.*Mixsigmin;
        
        if exist('OptSigs','Var')
            for ll=1:iNumofOptsigs
                    eval(sprintf('OptSig%dmin = rms_target./rms_min.*OptSig%dmin;',ll,ll));
                    eval(sprintf('OptSig%dmax = rms_target./rms_max.*OptSig%dmax;',ll,ll));
            end
        end
        % Use the modified SRMR and apply it to both EC processed signals
        [ratiomax(mm,kk), ~] = SRMRmod(Mixsigmax, fs_model, 'fast', 0, 'norm', 1, 'minCF', 4, 'maxCF', 32, 'single',0,'window',iBlockLen_EC);
        [ratiomin(mm,kk), ~] = SRMRmod(Mixsigmin, fs_model, 'fast', 0, 'norm', 1, 'minCF', 4, 'maxCF', 32, 'single',0,'window',iBlockLen_EC);
           
        if mm>2
             if mean([ratiomax(mm,kk) ratiomax(mm-1,kk)])>mean([ratiomin(mm,kk) ratiomin(mm-1,kk)]);
                Mixsigsyn = Mixsigmax; 
                if exist('OptSigs','Var')
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dsyn = OptSig%dmax;',ll,ll));
                    end
                end
                %fprintf('Max used in band %d\n',(kk)) KYLE O
            else
                Mixsigsyn = Mixsigmin;
                if exist('OptSigs','Var')
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dsyn = OptSig%dmin;',ll,ll));
                    end
                end
                %fprintf('Min used in band %d\n',(kk)) KYLE O
             end
         else
             
            if ratiomax(mm,kk)>ratiomin(mm,kk)
                Mixsigsyn = Mixsigmax; 
                if exist('OptSigs','Var')
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dsyn = OptSig%dmax;',ll,ll));
                    end
                end
                %fprintf('Max used in band %d\n',(kk)) KYLE O
            else
                Mixsigsyn = Mixsigmin;
                if exist('OptSigs','Var')
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dsyn = OptSig%dmin;',ll,ll));
                    end
                end
                %fprintf('Min used in band %d\n',(kk)) KYLE O
            end
        end
         
        % As the SII was used in the study of Hauth et al.(2020), the levels of the optional
        % signals (if present) are saved
        if exist('OptSigs','Var')
            for ll=1:iNumofOptsigs
                eval(sprintf('LevelOptSig%dmin(kk) = 20.*log10(rms(OptSig%dmin));',ll,ll));
                eval(sprintf('LevelOptSig%dmax(kk) = 20.*log10(rms(OptSig%dmax));',ll,ll));
                eval(sprintf('LevelOptSig%dsyn(kk) = 20.*log10(rms(OptSig%dsyn));',ll,ll));
                eval(sprintf('LevelOptSig%dL(kk) = 20.*log10(rms(OptSig%dL));',ll,ll));
                eval(sprintf('LevelOptSig%dR(kk) = 20.*log10(rms(OptSig%dR));',ll,ll));
            end
        end
         
        % Construct final signals for resynthesis
         if iNumofBlocks_EC==1
            Mixsigminfin(:,kk) = Mixsigminfin(:,kk) + Mixsigmin;
            Mixsigmaxfin(:,kk) = Mixsigmaxfin(:,kk) + Mixsigmax;
            Mixsigsynfin(:,kk) = Mixsigsynfin(:,kk)+ Mixsigsyn;
            MixsigLfin(:,kk)   = MixsigLfin(:,kk) + (rms_target./rmsL).* MixsigL;
            MixsigRfin(:,kk)   = MixsigRfin(:,kk) + (rms_target./rmsR).*MixsigR;
            
            if exist('OptSigs','Var')
                for ll=1:iNumofOptsigs
                    eval(sprintf('OptSig%dminfin(:,kk) = OptSig%dminfin(:,kk) + OptSig%dmin;',ll,ll,ll));
                    eval(sprintf('OptSig%dmaxfin(:,kk) = OptSig%dmaxfin(:,kk) + OptSig%dmax;',ll,ll,ll));
                    eval(sprintf('OptSig%dsynfin(:,kk) = OptSig%dsynfin(:,kk) + OptSig%dsyn;',ll,ll,ll));
                    eval(sprintf('OptSig%dLfin(:,kk) = OptSig%dLfin(:,kk) +(rms_target./rmsL).* OptSig%dL;',ll,ll,ll));
                    eval(sprintf('OptSig%dRfin(:,kk) = OptSig%dRfin(:,kk) +(rms_target./rmsR).* OptSig%dR;',ll,ll,ll));
                end
            end
     
         else
                     
            Mixsigminfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)...
            = Mixsigminfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+...
              Mixsigmin.*sqrt(window_EC);
        
            Mixsigmaxfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)...
            = Mixsigmaxfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+...
              Mixsigmax.*sqrt(window_EC);
          
            Mixsigsynfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)...
           = Mixsigsynfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+...
             Mixsigsyn.*sqrt(window_EC);
           
            MixsigLfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)...
            = MixsigLfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+...
              MixsigL.*sqrt(window_EC);
         
            MixsigRfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)...
            = MixsigRfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+...
              MixsigR.*sqrt(window_EC);
          
        % Do the same for optional signals if present
            if exist('OptSigs','Var')
                for ll=1:iNumofOptsigs
                    eval(sprintf('OptSig%dminfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk) = OptSig%dminfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+OptSig%dmin.*sqrt(window_EC);',ll,ll,ll));
                    eval(sprintf('OptSig%dmaxfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk) = OptSig%dmaxfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+OptSig%dmax.*sqrt(window_EC);',ll,ll,ll));
                    eval(sprintf('OptSig%dsynfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk) = OptSig%dsynfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+OptSig%dsyn.*sqrt(window_EC);',ll,ll,ll));
                    eval(sprintf('OptSig%dLfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk) = OptSig%dLfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+OptSig%dL.*sqrt(window_EC);',ll,ll,ll));
                    eval(sprintf('OptSig%dRfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk) = OptSig%dRfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+OptSig%dR.*sqrt(window_EC);',ll,ll,ll));
                end
            end
         end 
        end   
    else
        % Better ear processing:
        for mm = 1:iNumofBlocks_BE
            % Save left and right ear signal to temporary processing
            % variables.
            % if batch processing is used (whole signal is a single time
            % frame), save the whole signal.
            if iNumofBlocks_BE==1
                MixsigL = MixsigLbatch(:,kk);
                MixsigR = MixsigRbatch(:,kk);
                % Do the same for the optional signals.
                if exist('OptSigs','Var')
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dL = OptSig%dLbatch(:,kk);',ll,ll));
                        eval(sprintf('OptSig%dR = OptSig%dRbatch(:,kk);',ll,ll));
                    end
                end
            else
            % if short time processing is used, save only a signal portion.
                MixsigL = MixsigLbatch((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                          (mm-1)*iBlockshift_BE,kk).*sqrt(window_BE);
                MixsigR = MixsigRbatch((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                          (mm-1)*iBlockshift_BE,kk).*sqrt(window_BE);
                 % Do the same for the optional signals.
                if exist('OptSigs','Var')
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dL = OptSig%dLbatch((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk).*sqrt(window_BE);',ll,ll));
                        eval(sprintf('OptSig%dR = OptSig%dRbatch((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk).*sqrt(window_BE);',ll,ll));
                    end
                end  
            end
            
            rmsL = sqrt(mean(MixsigL.^2));
            rmsR = sqrt(mean(MixsigR.^2));
            rms_target = mean([rmsL' rmsR'],2);
            % Use the modified SRMR and apply it to the left and right ear
            % signal to identify the better ear.
            [ratioL, ~] = SRMRmod(rms_target./rmsL.*MixsigL, fs_model, 'fast', 0, 'norm', 1, 'minCF', 4, 'maxCF', 32, 'single',0,'window',iBlockLen_BE);
            [ratioR, ~] = SRMRmod(rms_target./rmsR.*MixsigR, fs_model, 'fast', 0, 'norm', 1, 'minCF', 4, 'maxCF', 32, 'single',0,'window',iBlockLen_BE);
           
            % Independent of the EC processing strategy, the better ear is
            % always selected, and the worse ear is always neglected. 
            % Therefore, the better ear channels are combined with all EC
            % outputs. (Min+BE, Max+BE, Syn+BE).
            if ratioL>=ratioR
                  % if the left ear channel has stronger modulation cues, 
                  % select it set it to the target level.
                  Mixsigmin  = rms_target./rmsL.*MixsigL;
                  Mixsigmax  = Mixsigmin;
                  Mixsigsyn = Mixsigmin;
                  % Do the same for optional signals
                  if exist('OptSigs','Var')
                    for ll=1:iNumofOptsigs  
                        eval(sprintf('OptSig%dsyn = rms_target./rmsL.*OptSig%dL;',ll,ll));
                        eval(sprintf('OptSig%dmin = rms_target./rmsL.*OptSig%dL;',ll,ll));
                        eval(sprintf('OptSig%dmax = rms_target./rmsL.*OptSig%dL;',ll,ll));
                    end
                  end
                  %disp('Left Ear used') KYLE O
            else
                  % if the right ear channel has stronger modulation cues, 
                  % select it and set it to the target level.
                  Mixsigmin  = rms_target./rmsR.*MixsigR;
                  Mixsigmax  = Mixsigmin;
                  Mixsigsyn = Mixsigmin;
                  % Do the same for optional signals
                  if exist('OptSigs','Var')
                    for ll=1:iNumofOptsigs 
                        eval(sprintf('OptSig%dsyn = rms_target./rmsR.*OptSig%dR;',ll,ll));
                        eval(sprintf('OptSig%dmin = rms_target./rmsR.*OptSig%dR;',ll,ll));
                        eval(sprintf('OptSig%dmax = rms_target./rmsR.*OptSig%dR;',ll,ll));
                    end
                  end
                  %disp('Right Ear used') KYLE O
            end
            % Because the SII was used in the study by Hauth et al.(2020),
            % the frequency specific levels of the optional signals are
            % calculated. They can later be used by the SII.
            % If no optional signals are present, the levels are not
            % calculated.
            if exist('OptSigs','Var')
                for ll=1:iNumofOptsigs
                    eval(sprintf('LevelOptSig%dmin(kk) = 20.*log10(rms(OptSig%dmin));',ll,ll));
                    eval(sprintf('LevelOptSig%dmax(kk) = 20.*log10(rms(OptSig%dmax));',ll,ll));
                    eval(sprintf('LevelOptSig%dsyn(kk) = 20.*log10(rms(OptSig%dsyn));',ll,ll));
                    eval(sprintf('LevelOptSig%dL(kk)   = 20.*log10(rms(OptSig%dL));',ll,ll));
                    eval(sprintf('LevelOptSig%dR(kk)   = 20.*log10(rms(OptSig%dR));',ll,ll));
                end
            end
             % save the channels in the final output matrices.
            if iNumofBlocks_BE==1
                % Minimized output
                Mixsigminfin(:,kk)=  Mixsigminfin(:,kk) + Mixsigmin;
                % Maximized output
                Mixsigmaxfin(:,kk)  = Mixsigmaxfin(:,kk)+ Mixsigmax;
                % Synthesized output
                Mixsigsynfin(:,kk) = Mixsigsynfin(:,kk) + Mixsigsyn;
                % Left ear output
                MixsigLfin(:,kk)  =  MixsigLfin(:,kk)   +(rms_target./rmsL).*MixsigL;
                % Right ear output
                MixsigRfin(:,kk)  =  MixsigRfin(:,kk)   +(rms_target./rmsR).*MixsigR;
                % Do the same for optional signals.
                if exist('OptSigs','Var')
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dminfin(:,kk)  = OptSig%dminfin(:,kk) + OptSig%dmin;',ll,ll,ll));
                        eval(sprintf('OptSig%dmaxfin(:,kk) = OptSig%dmaxfin(:,kk) + OptSig%dmax;',ll,ll,ll));
                        eval(sprintf('OptSig%dsynfin(:,kk) = OptSig%dsynfin(:,kk) + OptSig%dsyn;',ll,ll,ll));
                        eval(sprintf('OptSig%dLfin(:,kk)    = OptSig%dLfin(:,kk) +(rms_target./rmsL).* OptSig%dL;',ll,ll,ll));
                        eval(sprintf('OptSig%dRfin(:,kk)    = OptSig%dRfin(:,kk) +(rms_target./rmsR).* OptSig%dR;',ll,ll,ll));
                    end
                end
                % save channels to final output matrices ( if short time processing was used)
            else 
                %Mixed Signals
                % Mimizized
                Mixsigminfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)...
                = Mixsigminfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)+Mixsigmin.*sqrt(window_BE);
                % Maximimized
                Mixsigmaxfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)...
                = Mixsigmaxfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)+ Mixsigmax.*sqrt(window_BE);
                % Synthesized
                Mixsigsynfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)...
                = Mixsigsynfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)+ Mixsigsyn.*sqrt(window_BE);              
                 %Left
                MixsigLfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)...
                = MixsigLfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)+MixsigL.*sqrt(window_BE);
                %Right
                MixsigRfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)...
                = MixsigRfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)+MixsigR.*sqrt(window_BE);
           
                % Process Optional Signals accordingly:
                if exist('OptSigs','Var')
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dminfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk) = OptSig%dminfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk)+OptSig%dmin.*sqrt(window_BE);',ll,ll,ll));
                        eval(sprintf('OptSig%dmaxfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk) = OptSig%dmaxfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk)+OptSig%dmax.*sqrt(window_BE);',ll,ll,ll));
                        eval(sprintf('OptSig%dsynfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk) = OptSig%dsynfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk)+OptSig%dsyn.*sqrt(window_BE);',ll,ll,ll));
                        eval(sprintf('OptSig%dLfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk) = OptSig%dLfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk)+OptSig%dL.*sqrt(window_BE);',ll,ll,ll));
                        eval(sprintf('OptSig%dRfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk) = OptSig%dRfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk)+OptSig%dR.*sqrt(window_BE);',ll,ll,ll));
                    end
                end
            end
        end
    end
end
 %% Resynthesize time domain signals using gammatone filters and a delay of 10ms
synthesizer = gfb_synthesizer_new(analyzer1,0.01);%1/fs_model);%1/fs_model);
[Mixsigminfin, ~] = gfb_synthesizer_process(synthesizer, Mixsigminfin'+SomeImagSig);
[Mixsigmaxfin, ~] = gfb_synthesizer_process(synthesizer, Mixsigmaxfin');
[Mixsigsynfin, ~] = gfb_synthesizer_process(synthesizer, Mixsigsynfin');
[MixsigLfin, ~]   = gfb_synthesizer_process(synthesizer,MixsigLfin'+SomeImagSig);
[MixsigRfin, ~]   = gfb_synthesizer_process(synthesizer,MixsigRfin'+SomeImagSig);
% Do the same for optional signals
if exist('OptSigs','Var')
    for ll=1:iNumofOptsigs
        eval(sprintf('[OptSig%dminfin, ~] = gfb_synthesizer_process(synthesizer, transpose(OptSig%dminfin)+SomeImagSig);',ll,ll));
        eval(sprintf('[OptSig%dmaxfin, ~] = gfb_synthesizer_process(synthesizer, transpose(OptSig%dmaxfin));',ll,ll));
        eval(sprintf('[OptSig%dsynfin, ~] = gfb_synthesizer_process(synthesizer, transpose(OptSig%dsynfin));',ll,ll));
        eval(sprintf('[OptSig%dLfin, ~] = gfb_synthesizer_process(synthesizer, transpose(OptSig%dLfin)+SomeImagSig);',ll,ll));
        eval(sprintf('[OptSig%dRfin, ~] = gfb_synthesizer_process(synthesizer, transpose(OptSig%dRfin)+SomeImagSig);',ll,ll));
    end  
end

%% Resample signals if necessary
if fs_model~= fs
    Mixsigminfin  = resample(Mixsigminfin,fs,fs_model);
    Mixsigmaxfin  = resample(Mixsigmaxfin,fs,fs_model);
    Mixsigsynfin  = resample(Mixsigsynfin,fs,fs_model);
    MixsigLfin    = resample(MixsigLfin,fs,fs_model);
    MixsigRfin    = resample(MixsigRfin,fs,fs_model);
    
    if exist('OptSigs','Var')
        for ll=1:iNumofOptsigs
            eval(sprintf('[OptSig%dminfin, ~] = resample(OptSig%dminfin,fs,fs_model);',ll,ll));
            eval(sprintf('[OptSig%dmaxfin, ~] = resample(OptSig%dmaxfin,fs,fs_model);',ll,ll));
            eval(sprintf('[OptSig%dsynfin, ~] = resample(OptSig%dsynfin,fs,fs_model);',ll,ll));
            eval(sprintf('[OptSig%dLfin, ~]   = resample(OptSig%dLfin,fs,fs_model);',ll,ll));
            eval(sprintf('[OptSig%dRfin, ~]   = resample(OptSig%dRfinfs,fs_model);',ll,ll));
        end
    
    end
end
%% Save signals to output structure
out_struct.signals.MixsigSyn = Mixsigsynfin;
out_struct.signals.MixsigMin = Mixsigminfin;
out_struct.signals.MixsigMax = Mixsigmaxfin;
out_struct.signals.MixsigL = MixsigLfin;
out_struct.signals.MixsigR = MixsigRfin;
% Save also the optional signals to the output structure
    if exist('OptSigs','Var')
        for ll=1:iNumofOptsigs
            eval(sprintf('out_struct.signals.OptSig%dsyn = OptSig%dsynfin;',ll,ll));
            eval(sprintf('out_struct.levels.LevelOptSig%dsyn = LevelOptSig%dsyn;',ll,ll)); 
            eval(sprintf('out_struct.signals.OptSig%dmin = OptSig%dminfin;',ll,ll));
            eval(sprintf('out_struct.levels.LevelOptSig%dmin = LevelOptSig%dmin;',ll,ll));
            eval(sprintf('out_struct.signals.OptSig%dmax = OptSig%dmaxfin;',ll,ll));
            eval(sprintf('out_struct.levels.LevelOptSig%dmax = LevelOptSig%dmax;',ll,ll));
            eval(sprintf('out_struct.signals.OptSig%dL = OptSig%dLfin;',ll,ll));
            eval(sprintf('out_struct.levels.LevelOptSig%dL = LevelOptSig%dL;',ll,ll));
            eval(sprintf('out_struct.signals.OptSig%dR = OptSig%dRfin;',ll,ll));
            eval(sprintf('out_struct.levels.LevelOptSig%dR = LevelOptSig%dR;',ll,ll));
        end
    end
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