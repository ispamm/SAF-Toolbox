function [F, y, e] = AF_LMS_WPOLY_F(F, x, d)

% This function implements the adaptation of an FIR Wiener polynomial 
% structure by using the LMS adaptive algorithm.
%
% Details can be found in:
% 
%
%
% USAGE:
%  [F, y, e] =  AF_LMS_WPOLY_F(F, x, d)
%
% Input Arguments:
%   F       : AF structure
%   x       : array of input signal
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
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------


if nargin==0, help AF_LMS; return; end

M = F.M;                                     % Length of the linear filter

% Linear Filter -----------------------------------------------------------
F.xw(2:M)        = F.xw(1:M-1);              % Shift the input delay-line 
F.xw(1)          = x;                        % Load a new input into the delay-line
s                = F.w'*F.xw;                % Linear combiner output
[y, F.af]        = ActFunc(s,F.af);          % Nonlinear output
e                = d - y;                    % Error evaluation
ee               = e*dActFunc(s, F.af);      % Error multiplied by the derivative
F.w              = F.w + F.mu*conj(ee)*F.xw; % Updating the filter taps

% Polynomial Nonlinearity -------------------------------------------------
if F.af.aftype == 3    % 3 = Polynomial
    e_av = F.mQ*e;
    F.af.Q = F.af.Q + e_av*F.af.g';          % Updating the polynomial coefficients
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------