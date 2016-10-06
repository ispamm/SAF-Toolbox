function  [F, y, s] =  FW_HSPL_F(F, x)

% This function evaluate the output of a Hammerstein nonlinear structure 
% implemented by a HSAF architecture.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Hammerstein 
% Uniform Cubic Spline Adaptive Filters: Learning and Convergence 
% Properties", Signal Processing, Vol. 100, pp. 112-123, July 2014.
%
% USAGE:
%   [F, y(n), s] = FW_HSPL_F(F, x(n))
%
% Input Arguments:
%   F       : ADAPTIVE FILTER STRUCT
%   x       : x[n] input signal sample
%
% Output Arguments:
%   F       : ADAPTIVE FILTER STRUCT
%   y       : y[n] output signal sample
%   s       : s linear combiner input array
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2013/05/15$
% $Revision: 1.1$  $Date: 2016/07/31$
% -------------------------------------------------------------------------
% License to use and modify this code is granted freely without warranty to
% all, as long as the original authors are referenced and attributed as such.
% The original authors maintain the right to be solely associated with this
% work.
% -------------------------------------------------------------------------
% © 2013
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------


M = F.M;                      % Length of the linear filter

F.xw(2 : M) = F.xw(1 : M-1);  % Shift of the input delay-line
F.xw(1) = x;                  % Load a new input into the delay-line 

s = zeros(M,1);               % Linear combiner input array
for i=1:M,
    [s(i),F.af] = ActFunc(F.xw(i),F.af); % Evaluating the nonlinearity output
end
y = s'*F.w;                   % Filter output

% -------------------------------------------------------------------------
% © 2013
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------