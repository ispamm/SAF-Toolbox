function [F, y, s, r] = FW_Sandwich1SPL_F(F, x)

% This function evaluate the output of a Sandwich 1 nonlinear structure 
% implemented by a S1SAF architecture.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Novel Cascade 
% Spline Architectures for the Identification of Nonlinear Systems", 
% IEEE Transactions on Circuits and Systems---I: Regular Papers, Vol. 62,
% No. 7, pp. 1825-1835, July 2015.
%
% USAGE:
%   [F, y(n), s(n), r(n)] = FW_Sandwich1SPL_F(F, x(n))
%
% Input Arguments:
%   F       : ADAPTIVE FILTER STRUCT
%   x       : x[n] input signal sample
%
% Output Arguments:
%   F       : ADAPTIVE FILTER STRUCT
%   y       : y[n] output signal sample
%   s       : s[n] first nonlinear output
%   r       : r[n] linear filter output
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2014/10/15$
% $Revision: 1.1$  $Date: 2016/08/18$
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


M = F.M;                       % Length of the linear filter

F.xw(2 : M) = F.xw(1 : M-1);   % Shift the input filter delay-line
F.xw(1) = x;                   % Load a new input into the delay-line 

s = zeros(M,1);
for i=1:M,
    [s(i),F.af1] = ActFunc(F.xw(i),F.af1);  % Output of first nonlinearity in Eq. (11)
end
r = s'*F.w;                    % Linear filter output in Eq. (12)
[y,F.af2] = ActFunc(r,F.af2);  % System output in Eq. (13)

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------