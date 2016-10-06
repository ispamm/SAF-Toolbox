function [F, y, s, r] = FW_Sandwich2SPL_F(F, x)

% This function evaluate the output of a Sandwich 2 nonlinear structure 
% implemented by a S2SAF architecture.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Novel Cascade 
% Spline Architectures for the Identification of Nonlinear Systems", 
% IEEE Transactions on Circuits and Systems---I: Regular Papers, Vol. 62,
% No. 7, pp. 1825-1835, July 2015.
%
% USAGE:
%   [F, y(n), s(n), r(n)] = FW_Sandwich2SPL_F(F, x(n))
%
% Input Arguments:
%   F       : ADAPTIVE FILTER STRUCT
%   x       : x[n] input signal sample
%
% Output Arguments:
%   F       : ADAPTIVE FILTER STRUCT
%   y       : y[n] output signal sample
%   s       : s[n] first linear filter output
%   r       : r[n] nonlinearity output
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


M1 = F.M1;                        % Length of the first linear filter
M2 = F.M2;                        % Length of the second linear filter

F.xw1(2 : M1) = F.xw1(1 : M1-1);  % Shift the input filter 1 delay-line
F.xw1(1) = x;                     % Load a new input into the delay-line 
s = F.xw1'*F.w1;                  % First linear filter output
[r,F.af] = ActFunc(s,F.af);       % Nonlinearity output
F.xw2(2 : M2) = F.xw2(1 : M2-1);  % Shift the input filter 2 delay-line
F.xw2(1) = r;                     % Load a new input into the delay-line
y = F.xw2'*F.w2;                  % System output
  
% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------