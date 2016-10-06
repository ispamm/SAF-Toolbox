function lms_af = create_Sandwich1SPL_lms_adaptive_filter_1(M, mu, mQ1, mQ2, delta, af1, af2)

% This function creates and initializes the structure implementing a S1SAF 
% architecture and adapted by the LMS adaptive algorithm.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Novel Cascade 
% Spline Architectures for the Identification of Nonlinear Systems", 
% IEEE Transactions on Circuits and Systems---I: Regular Papers, Vol. 62,
% No. 7, pp. 1825-1835, July 2015.
%
% USAGE:
% lms_af =  create_Sandwich1SPL_lms_adaptive_filter_1(M, mu, mQ1, muQ2, delta, af1, af2)
%   - M: length of the linear filter
%   - mu: step-size of the lineaar filter
%   - mQ1: step-size of the first nonlinearity
%   - mQ2: step-size of the second nonlinearity
%   - delta: regularization parameter
%   - af1: first spline nonlinearity structure
%   - af2: second spline nonlinearity structure
%   - lms_af: structure of the created SAF filter
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


% LMS Sandwich 1 adaptive filter definition and initialization ------------
w = zeros(M,1);     % Linear filter taps  
w(floor(M/2)) = 1;  % Adaptive filter weights I.C.
dI = delta;         % Regularizing parameter
xd = zeros(M,1);    % Buffer of the desired signal
xw = zeros(M,1);    % Buffer of the filter status
% -------------------------------------------------------------------------
lms_af = struct('M',M, 'mu',mu,'mQ1',mQ1, 'mQ2',mQ2, 'dI',dI, 'w',w, ...
    'xd',xd, 'xw',xw, 'af1',af1, 'af2',af2);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------