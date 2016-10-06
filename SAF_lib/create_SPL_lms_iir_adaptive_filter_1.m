function  lms_af = create_SPL_lms_iir_adaptive_filter_1(M, N, mu, mQ, delta, af)

% This function creates and initializes the structure implementing a IIR 
% WSAF architecture and adapted by the LMS adaptive algorithm.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear System
% Identification using IIR Spline Adaptive Filters", Signal Processing, 
% Vol. 108, pp. 30-35, March 2015.
%
% USAGE:
% lms_af =  create_SPL_lms_iir_adaptive_filter_1(M, N, mu, mQ, delta, af)
%   - M: length of the MA part of the linear filter
%   - N: length of the AR part of the linear filter
%   - mu: step-size of the lineaar filter
%   - mQ: step-size of the nonlinearity
%   - delta: regularization parameter
%   - af: spline nonlinearity structure
%   - lms_af: structure of the created SAF filter
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


% LMS IIR SAF adaptive filter definition and initialization ---------------
Mt    = M + N;          % Total number of parametrs
b     = zeros(M,1);     % MA taps
a     = zeros(N,1);     % AR taps
w     = zeros(Mt,1);    % Complete AF taps (ARMA model)  
beta  = zeros(M,N+1);   % Derivatives with respect MA coefficients
alpha = zeros(N,N+1);   % Derivatives with respect AR coefficients
dI = delta;             % Regularizing parameter
xd = zeros(Mt,1);       % Buffer of the desired signal
xw = zeros(M,1);        % Buffer of the filter status
sw = zeros(N+1,1);      % Buffer of the IIR status
% -------------------------------------------------------------------------
lms_af = struct('M',M, 'N',N, 'mu',mu,'mQ',mQ, 'dI',dI, 'w',w, 'b',b, 'a',a, ...
    'beta',beta, 'alpha',alpha, 'xd',xd, 'xw',xw, 'sw',sw, 'af',af);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------