function [F, y, e] = AF_APA_WSPL_IIR_F(F, x, d)

% This function implements the adaptation of an IIR Wiener SAF (WSAF) 
% structure by using the APA adaptive algorithm.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear System
% Identification using IIR Spline Adaptive Filters", Signal Processing, 
% Vol. 108, pp. 30-35, March 2015.
%
% USAGE:
%  [F, y, e] =  AF_APA_WSPL_IIR_F(F, x, d)
%
% Input Arguments:
%   F       : IIR WSAF structurs
%   x       : array of input signal
%   d       : array of desired signal
%  
% Output Arguments:
%   F       : IIR WSAF structure
%   y       : array of output signal   
%   e       : array of error signal    
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2014/04/10$
% $Revision: 1.1$  $Date: 2016/07/31$
% -------------------------------------------------------------------------
% License to use and modify this code is granted freely without warranty to
% all, as long as the original authors are referenced and attributed as such.
% The original authors maintain the right to be solely associated with this
% work.
% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------

% !!! WORK IN PROGRESS !!!

M = F.M;   % Length of the MA part
N = F.N;   % Length of the AR part
K = F.K;   % APA order
yy = zeros(K,1);

% Covar Matrix Up Date ---------------------------------------------------- 
F.xw(2 : M) = F.xw(1 : M-1);     % Shift of the input delay-line
F.xw(1) = x;                     % Load a new input into the delay-line  
for i = 1:K-1    
    F.X(i,1:M) = F.X(i+1,1:M);   % Shift the covarince matrix X 
end
F.X(K,1:M) = F.xw(1:M);          % Fill the covariance matrix X 
    
% Desired-output buffer for APA -------------------------------------------
F.xd(1:K-1) = F.xd(2:K);        % Shift of the input delay-line  
F.xd(K) = d ;                   % Load a new desired output into the delay-line 

ss = F.X*F.w;                   % (K,1) = (K,M) x (M,1) output array 
for k=1:K
    [yy(k),F.af] = ActFunc(ss(k),F.af);
end 
   
ee = F.xd - yy;    % (K,1) = (K,1) - (K,1) error array 
e = ee(K);         % A priori output error sample
y = yy(K);         % Output sample
    
% Learning ------------------------------------------------------------ 
for k=1:K  %  err*phi'(s) 
    ee(k) = ee(k)*dActFunc(y,F.af);
end 
  
% APA weights updating ----------------------------------------------------
F.w = F.w + F.mu*(F.X'/(F.dI + F.X*F.X'))*ee;  % (M,1) = (M,1) + (M,K) x (K,1)  

% LMS LUT control points updating -----------------------------------------
if F.af.aftype > 1      % Kind act. f. -1 0 1 2 4 5 6 7 ... 20  */
    e_av = (F.mQ*e)/1;
    ii = F.af.uIndex :  F.af.uIndex + F.af.P;   % P = Spline order
    F.af.Q(ii) = F.af.Q(ii) + e_av*F.af.g';
end
    
% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------