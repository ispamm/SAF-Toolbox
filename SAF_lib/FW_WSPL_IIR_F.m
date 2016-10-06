function  [F, y] =  FW_WSPL_IIR_F(F, x)

% This function evaluate the output of a IIR Wiener nonlinear structure 
% implemented by a IIR WSAF architecture.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear System
% Identification using IIR Spline Adaptive Filters", Signal Processing, 
% Vol. 108, pp. 30-35, March 2015.
%
% USAGE:
%   [F, y(n)] = FW_WSPL_IIR_F(F, x(n))
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


M = F.M;                                 % Length of the MA part
N = F.N;                                 % Length of the AR part

F.xw(2 : M) = F.xw(1 : M-1);             % Shift the input delay-line
F.xw(1) = x;                             % Load a new input into the delay-line 
F.sw(2:N+1) = F.sw(1:N);                 % Shift the internal delay-line
F.sw(1) = F.xw'*F.b + F.sw(2:N+1)'*F.a;  % Load a new input into the internal delay-line with the output of an IIR filter
[y,F.af] = ActFunc(F.sw(1),F.af);        % Output of the IIR WSAF

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------