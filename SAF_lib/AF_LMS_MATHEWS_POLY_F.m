function [F, y, e] = AF_LMS_MATHEWS_POLY_F(F, x, d)

% This function implements the adaptive LNL sandwich model proposed in:
% J. Jeraj and V. J. Mathews, "A stable adaptive Hammerstein filter employing
% partial orthogonalization of the input signals", IEEE Transactions on
% Signal Processing, Vol. 54, N. 4, pp. 1412-1420, April 2006.
%
% It has been used for comparisons in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Novel Cascade 
% Spline Architectures for the Identification of Nonlinear Systems", 
% IEEE Transactions on Circuits and Systems---I: Regular Papers, Vol. 62,
% No. 7, pp. 1825-1835, July 2015.
%
% USAGE:
%  [F, y, e] =  AF_LMS_MATHEWS_POLY_F(F, x, d)
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
% $Revision: 1.0$  $Date: 2014/10/15$
% $Revision: 1.1$  $Date: 2016/08/31$
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


M = F.M;                             % Length of the linear filter
P = F.af.Pord;                       % Polynomial order

theta = [F.w(2:M).' F.af.Q.'].';
Lambda = diag(F.mu*ones(1,M+P-1));

F.xw(2 : M) = F.xw(1 : M-1);         % Shift the filter delay-line
[F.xw(1), F.af] = ActFunc(x, F.af);  % Load a new input into the delay-line  
F.af.gM(:,2:M) = F.af.gM(:,1:M-1);   % Shift the filter delay-line
F.af.gM(:,1) = F.af.g.';             % Load a new sample into the delay-line

H = [F.xw(2:M).' F.af.g].';             % Construct the matrix H
Psi = [F.xw(2:M).' (F.af.gM*F.w).'].';  % Construct the matrix Psi

y = H.'*theta;                       % Output evaluation
e = d - y;                           % Error evaluation

% LMS weights and Polynomial coefficients update --------------------------
theta = theta + Lambda*Psi*inv(F.dI + H.'*Lambda*Psi)*e;  % Update all parameters

F.w(2:M) = theta(1:M-1);             % Linear filter weights
F.w(1) = 1;                          % First tap is always 1
F.af.Q = theta(M:M+P-1);             % Polynomials coefficients
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------