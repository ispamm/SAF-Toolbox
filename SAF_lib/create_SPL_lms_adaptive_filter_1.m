function lms_af = create_SPL_lms_adaptive_filter_1(M, mu, mQ, delta, af)

% This function creates and initializes the structure implementing a WSAF 
% architecture and adapted by the LMS adaptive algorithm.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear Spline 
% Adaptive Filtering", Signal Processing, Vol. 93, No. 4, pp. 772-783, 
% April 2013.
%
% USAGE:
% lms_af =  create_SPL_lms_adaptive_filter_1(M, mu, mQ, delta, af)
%   - M: length of the linear filter
%   - mu: step-size of the lineaar filter
%   - mQ: step-size of the nonlinearity
%   - delta: regularization parameter
%   - af: spline nonlinearity structure
%   - lms_af: structure of the created SAF filter
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


% LMS SAF adaptive filter definition and initialization -------------------
w = zeros(M,1);     % Linear filter taps  
w(floor(M/2)) = 1;  % Adaptive filter weights i.c.
dI = delta;         % Regularizing parameter
xd = zeros(M,1);    % Buffer of the desired signal
xw = zeros(M,1);    % Buffer of the filter status
% ---------------------------------------------------------------------------------
lms_af = struct('M',M, 'mu',mu,'mQ',mQ, 'dI',dI, 'w',w, 'xd',xd, 'xw',xw, 'af',af);
% ---------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------