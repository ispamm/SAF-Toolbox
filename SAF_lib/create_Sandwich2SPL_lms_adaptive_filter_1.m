function lms_af = create_Sandwich2SPL_lms_adaptive_filter_1(M1, M2, mu1, mu2, mQ, delta, af)

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
% lms_af =  create_Sandwich2SPL_lms_adaptive_filter_1(M, mu, mQ1, muQ2, delta, af1, af2)
%   - M1: length of the first linear filter
%   - M2: length of the second linear filter
%   - mu1: step-size of the first lineaar filter
%   - mu2: step-size of the second lineaar filter
%   - mQ: step-size of the nonlinearity
%   - delta: regularization parameter
%   - af: spline nonlinearity structure
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

% LMS Sandwich 2 adaptive filter definition and initialization ------------
w1  = zeros(M1,1);    % Linear filter 1 taps  
w2  = zeros(M2,1);    % Linear filter 2 taps  
w1(floor(M1/2)) = 1;  % Adaptive filter weights i.c.
w2(floor(M2/2)) = 1;  % Adaptive filter weights i.c.
dI  = delta;          % Regularizing parameter
xd  = zeros(M2,1);    % Buffer of the desired signal
xw1 = zeros(M1,1);    % Buffer of the filter 1 status
xw2 = zeros(M2,1);    % Buffer of the filter 2 status
xs  = zeros(M2,1);    % Buffer of the nonlinearity
X   = zeros(M1,M2);   % Buffer matrix of the past inputs delay lines
rp  = zeros(M2,1);    % Buffer of the past nonlinearity derivatives
% -------------------------------------------------------------------------
lms_af = struct('M1',M1, 'M2',M2, 'mu1',mu1, 'mu2',mu2, 'mQ',mQ, 'dI',dI, ...
    'w1',w1, 'w2',w2, 'xd',xd, 'xw1',xw1, 'xw2',xw2, 'X',X, 'rp',rp, 'xs',xs, 'af',af);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------