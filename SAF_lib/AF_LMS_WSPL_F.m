function [F, y, e] = AF_LMS_WSPL_F(F, x, d)

% This function implements the adaptation of a Wiener SAF (WSAF) structure 
% by using the LMS adaptive algorithm.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear Spline 
% Adaptive Filtering", Signal Processing, Vol. 93, No. 4, pp. 772-783, 
% April 2013.
%
% USAGE:
%  [F, y, e] =  AF_LMS_WSPL_F(F, x, d)
%
% Input Arguments:
%   F       : WSAF structurs
%   x       : array of input signal
%   d       : array of desired signal
%  
% Output Arguments:
%   F       : WSAF structure
%   y       : array of output signal   
%   e       : array of error signal    
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2012/01/22$
% $Revision: 1.1$  $Date: 2016/07/31$
% -------------------------------------------------------------------------
% License to use and modify this code is granted freely without warranty to
% all, as long as the original authors are referenced and attributed as such.
% The original authors maintain the right to be solely associated with this
% work.
% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------

M         = F.M;                % Length of the linear filter
F.xw(2:M) = F.xw(1:M-1);        % Shift of the filter delay-line 
F.xw(1)   = x;                  % Load a new input into the delay-line
s         = F.w'*F.xw;          % Otput of the linear filter
[y,F.af]  = ActFunc(s,F.af);    % Output of the nonlinearity
e         = d - y;              % Error evaluation
ee        = e*dActFunc(s,F.af); % Error multiplied by the derivative

% LMS weights and control points update -----------------------------------
F.w      = F.w + F.mu*conj(ee)*F.xw;   % cpx LMS in Eq. (15)

if F.af.aftype > 1   % Kind act. f. -1 0 1 2 3 4 5 6 7 ... 20  */
    e_av = F.mQ*e;
    ii = F.af.uIndex : F.af.uIndex + F.af.P;   % P = Spline order
    F.af.Q(ii) = F.af.Q(ii) + e_av*F.af.g';    % LMS in Eq. (16)
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------