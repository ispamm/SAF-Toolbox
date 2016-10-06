function  [F, y, s] =  FW_WPOLY_F(F, x)

% This function evaluate the output of a FIR Wiener polynomial nonlinear 
% structure implemented by a polynomial AF.
%
% Details can be found in:
% 
%
%
% USAGE:
%   [F, y(n), s(n)] = FW_WPOLY_F(F, x(n))
%
% Input Arguments:
%   F       : ADAPTIVE FILTER STRUCT
%   x       : x[n] input signal sample
%   s       : s[n] linear combiner output
%
% Output Arguments:
%   F       : ADAPTIVE FILTER STRUCT
%   y       : y[n] output signal sample
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


M = F.M;                      % Length of the linear filter

F.xw(2 : M) = F.xw(1 : M-1);  % Shift the input delay-line
F.xw(1) = x;                  % Load a new input into the delay-line 
   
s = F.xw.'*F.w;               % Linear filter output
y = ActFunc(s,F.af);          % Nonlinear output
 
% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------