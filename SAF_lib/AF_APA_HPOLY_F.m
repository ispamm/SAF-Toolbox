function  [F, y, e] =  AF_APA_HPOLY_F(F, x, d)

% This function implements the adaptation of an FIR Hammerstein polynomial 
% structure by using the APA adaptive algorithm.
%
% Details can be found in:
% 
%
%
% USAGE:
%  [F, y, e] =  AF_APA_HPOLY_F(F, x, d)
%
% Input Arguments:
%   F       : AF structure
%   u       : array of input signal
%   d       : array of desired signal
%  
% Output Arguments:
%   F       : AF structure
%   y       : array of output signal   
%   e       : array of error signal    
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2014/04/14$
% $Revision: 1.1$  $Date: 2016/08/01$
% -------------------------------------------------------------------------
% License to use and modify this code is granted freely without warranty to
% all, as long as the original authors are referenced and attributed as such.
% The original authors maintain the right to be solely associated with this
% work.
% -------------------------------------------------------------------------
% � 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------


M = F.M;              % Length of the linear filter
K = F.K;              % APA order

% Covariance Matrix Update ------------------------------------------------ 
F.xw(2 : M) = F.xw(1 : M-1);         % Shift the input delay-line
[F.xw(1), F.af] = ActFunc(x, F.af);  % Load a new input into the delay-line  
for i = 1:K-1                        % Shift the covarince matrix X
    F.X(i,1:M) = F.X(i+1,1:M); 
end
F.X(K,1:M) = F.xw(1:M);              % Fill the covariance matrix X 

% Desired-output buffer for APA -------------------------------------------
F.xd(1:K-1) = F.xd(2:K);             % Shift desired delay-line  
F.xd(K) = d ;                        % Load a new desired sample into the delay-line 

yy = F.X*F.w;                        % (K,1) = (K,M) x (M,1) output array 
ee = F.xd - yy;                      % (K,1) = (K,1) - (K,1) error array 
e = ee(K);                           % A priori output error sample
y = yy(K);                           % Output sample

% Linear Filter -----------------------------------------------------------
F.w = F.w + F.mu*(F.X'/(F.dI + F.X*F.X'))*ee;  % Updating the filter taps 
F.w(1) = 1;                                    % If first tap should be 1

% Polynomial Nonlinearity -------------------------------------------------
if (F.af.aftype == 3)        % 3 = Polynomial
    F.af.gM(:,2:M) = F.af.gM(:,1:M-1);         % Shift the U matrix buffer
    F.af.gM(:,1) = F.af.g.';                   % Load a new vector in matrix U
    u = F.af.gM*F.w;                           % Evaluating u

    e_av = F.mQ*e;     
    F.af.Q = F.af.Q + e_av*u;                  % Updating the polynomial coefficients
end 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% � 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------