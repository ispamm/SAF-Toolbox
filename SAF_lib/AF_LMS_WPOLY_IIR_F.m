function   [F, y, e] = AF_LMS_WPOLY_IIR_F(F, x, d)

% This function implements the adaptation of an IIR Polynomial Wiener (WPOLY) 
% structure by using the LMS adaptive algorithm.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear System
% Identification using IIR Spline Adaptive Filters", Signal Processing, 
% Vol. 108, pp. 30-35, March 2015.
%
% USAGE:
%  [F, y, e] =  AF_LMS_WPOLY_IIR_F(F, x, d)
%
% Input Arguments:
%   F       : IIR WPOLY structurs
%   x       : array of input signal
%   d       : array of desired signal
%  
% Output Arguments:
%   F       : IIR WPOLY structure
%   y       : array of output signal   
%   e       : array of error signal    
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2014/04/14$
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


if nargin==0, help AF_LMS; return; end

M = F.M;  % Length of the MA part
N = F.N;  % Length of the AR part

% Linear Filter -----------------------------------------------------------
F.xw(2:M)        = F.xw(1:M-1);                        % Shift the input delay-line
F.xw(1)          = x;                                  % Load a new input into the delay-line
F.sw(2:N+1)      = F.sw(1:N);                          % Shift the internal delay-line
F.sw(1)          = F.b'*F.xw + F.a'*F.sw(2:N+1);       % Load a new input into the internal delay-line with the output of an IIR filter
[y, F.af]        = ActFunc(F.sw(1),F.af);              % Output of the IIR Polynomial AF
e                = d - y;                              % Error evaluation
ee               = e*dActFunc(y, F.af);                % Error multiplied by the nonlinear derivative
F.w              = [F.b.'  F.a.'].';                   % ARMA parameters
F.beta(:,2:N+1)  = F.beta(:,1:N);                      % Shift the beta delay-line
F.beta(:,1)      = F.xw + F.beta(:,2:N+1)*F.a;         % Updating the beta delay-line
F.alpha(:,2:N+1) = F.alpha(:,1:N);                     % Shift the alpha delay-line
F.alpha(:,1)     = F.sw(2:N+1) + F.alpha(:,2:N+1)*F.a; % Updating the alpha delay-line
eta              = [F.beta(:,1).'  F.alpha(:,1).'].';  % Constructing the eta vector
F.w              = F.w + F.mu*conj(ee)*eta;            % Updating all the parameters
F.b              = F.w(1:M);                           % Updating the MA parameters
F.a              = F.w(M+1:M+N);                       % Updating the AR parameters

% Polynomial Nonlinearity -------------------------------------------------
if F.af.aftype == 3      % 3 = Polynomial
    e_av = F.mQ*e;
    F.af.Q = F.af.Q + e_av*F.af.g';   % Coefficients update
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------