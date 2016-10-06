function  [F, y] =  FW_WSPL_F(F, x)

% This function evaluate the output of a Wiener nonlinear structure 
% implemented by a WSAF architecture.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear Spline 
% Adaptive Filtering", Signal Processing, Vol. 93, No. 4, pp. 772-783, 
% April 2013.
%
% USAGE:
%   [F, y(n)] = FW_WSPL_F(F, x(n))
%
% Input Arguments:
%   F       : ADAPTIVE FILTER STRUCT
%   x       : x[n] input signal sample
%
% Output Arguments:
%   F       : ADAPTIVE FILTER STRUCT
%   y       : y[n] output signal sample
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


M = F.M;     % Length of the linear filter

F.xw(2 : M) = F.xw(1 : M-1);   % Shift the input delay-line
F.xw(1) = x;                   % Load a new input into the delay-line 

[y, F.af] = ActFunc(F.xw'*F.w, F.af); % Output of the nonlinearity, Eq. (8)


% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------