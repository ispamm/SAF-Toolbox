function   [F, y, e] = AF_LMS_WSPL_IIR_F(F, x, d)

% This function implements the adaptation of an IIR Wiener SAF (WSAF) 
% structure by using the LMS adaptive algorithm.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear System
% Identification using IIR Spline Adaptive Filters", Signal Processing, 
% Vol. 108, pp. 30-35, March 2015.
%
% USAGE:
%  [F, y, e] =  AF_LMS_WSPL_IIR_F(F, x, d)
%
% Input Arguments:
%   F       : IIR WSAF structure
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


if nargin==0, help AF_LMS; return; end

M = F.M;  % Length of the MA part
N = F.N;  % Length of the AR part

% Linear Filter -----------------------------------------------------------
F.xw(2:M)        = F.xw(1:M-1);                   % Shift the input delay-line
F.xw(1)          = x;                             % Load a new input into the delay-line 
F.sw(2:N+1)      = F.sw(1:N);                     % Shift the internal delay-line
F.sw(1)          = F.b'*F.xw + F.a'*F.sw(2:N+1);  % Load a new input into the internal delay-line with the output of an IIR filter
[y, F.af]        = ActFunc(F.sw(1),F.af);         % Output of the IIR WSAF
e                = d - y;                         % Error evaluation
ee               = e*dActFunc(y, F.af);           % Error multiplied by the nonlinear derivative
F.w              = [F.b.'  F.a.'].';              % ARMA parameters
F.beta(:,2:N+1)  = F.beta(:,1:N);                 % Shift the beta delay-line
F.beta(:,1)      = F.xw + F.beta(:,2:N+1)*F.a;    % Updating the beta delay-line in Eq. (14)
F.alpha(:,2:N+1) = F.alpha(:,1:N);                % Shift the alpha delay-line
F.alpha(:,1)     = F.sw(2:N+1) + F.alpha(:,2:N+1)*F.a; % Updating the alpha delay-line in Eq. (15)
eta              = [F.beta(:,1).'  F.alpha(:,1).'].';  % Constructing the eta vector in Eq. (16)
F.w              = F.w + F.mu*conj(ee)*eta;       % Updating all the parameters as in Eq. (17)
F.b              = F.w(1:M);                      % Updating the MA parameters
F.a              = F.w(M+1:M+N);                  % Updating the AR parameters

% Spline Nonlinearity -----------------------------------------------------
if F.af.aftype > 1    % Kind act. f. -1 0 1 2 3 4 5 6 7 ... 20  */
    e_av = F.mQ*e;
    ii = F.af.uIndex : F.af.uIndex + F.af.P;      % P = Spline order
    for  m = 1:1
        F.af.Q(ii) = F.af.Q(ii) + e_av*F.af.g';   % LMS in Eq. (19)
    end
end
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------