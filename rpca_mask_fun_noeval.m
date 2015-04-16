function [wavoutVoc,wavoutAcc]=rpca_mask_fun_noeval(wavinmix,param)
    %% parameters
    if nargin>1
        lambda = param.lambda;
        nFFT = param.nFFT;
        winsize = param.windowsize;
        masktype = param.masktype;
        gain = param.gain;
        power = param.power;
        % Fs= parm.fs;
        % outputname = parm.outname;
    else
        lambda=1;
        nFFT=1024;
        winsize=1024;
        masktype=1; %1: binary mask, 2: no mask
        gain=1;
        power=1;
    end

    hop = winsize/4;
    scf = 2/3;
    S = scf * stft(wavinmix, nFFT ,winsize, hop);

   %% use inexact_alm_rpca to run RPCA
   oldpath=addpath(genpath('inexact_alm_rpca'));
    try                
        [A_mag, E_mag] = inexact_alm_rpca(abs(S).^power',lambda/sqrt(max(size(S))));
        PHASE = angle(S');            
    catch 
        [A_mag, E_mag] = inexact_alm_rpca(abs(S).^power,lambda/sqrt(max(size(S))));
        PHASE = angle(S);
    end
    
    A = A_mag.*exp(1i.*PHASE);
    E = E_mag.*exp(1i.*PHASE);
    path(oldpath);
    %% binary mask, no mask
    switch masktype                         
      case 1 % binary mask + median filter
        m= double(abs(E)> (gain*abs(A)));                  
        try  
            Emask =m.*S;
            Amask= S-Emask;
        catch
            Emask =m.*S';
            Amask= S'-Emask;
        end        
      case 2 % no mask
        Emask=E;
        Amask=A;
      otherwise 
          fprintf('masktype error\n');
    end

    %% do istft
    try 
        wavoutVoc = istft(Emask', nFFT ,winsize, hop)';   
        wavoutAcc = istft(Amask', nFFT ,winsize, hop)';
    catch 
        wavoutVoc = istft(Emask, nFFT ,winsize, hop)';   
        wavoutAcc = istft(Amask, nFFT ,winsize, hop)';
    end

    wavoutVoc=wavoutVoc/max(abs(wavoutVoc));
%     wavwrite(wavoutE,Fs,[outputname,'_E']);

    wavoutAcc=wavoutAcc/max(abs(wavoutAcc));
%     wavwrite(wavoutA,Fs,[outputname,'_A']);

    %% evaluate
%     if length(wavoutA)==length(wavinA)
% 
%         sep = [wavoutA , wavoutE]';
%         orig = [wavinA , wavinE]';
% 
%         for i = 1:size( sep, 1)
%                [e1,e2,e3] = bss_decomp_gain( sep(i,:), i, orig);
%                [sdr(i),sir(i),sar(i)] = bss_crit( e1, e2, e3);
%         end
%     else
%         minlength=min( length(wavoutE), length(wavinE) );
% 
%         sep = [wavoutA(1:minlength) , wavoutE(1:minlength)]';
%         orig = [wavinA(1:minlength) , wavinE(1:minlength)]';
% 
%         for i = 1:size( sep, 1)
%                [e1,e2,e3] = bss_decomp_gain( sep(i,:), i, orig);
%                [sdr(i),sir(i),sar(i)] = bss_crit( e1, e2, e3);
%         end
%     end
% 
%     Parms.SDR=sdr(2);
%     Parms.SIR=sir(2);
%     Parms.SAR=sar(2);
