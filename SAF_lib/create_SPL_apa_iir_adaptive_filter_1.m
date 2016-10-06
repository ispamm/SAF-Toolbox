function apa_af = create_SPL_apa_iir_adaptive_filter_1(M, N, K, mu, mQ, delta, af)

% This function creates and initializes the structure implementing a IIR 
% WSAF architecture and adapted by the APA adaptive algorithm.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear System
% Identification using IIR Spline Adaptive Filters", Signal Processing, 
% Vol. 108, pp. 30-35, March 2015.
%
% USAGE:
% apa_af = create_SPL_apa_iir_adaptive_filter_1(M, N, mu, mQ, delta, af)
%   - M: length of the MA part of the linear filter
%   - N: length of the AR part of the linear filter
%   - K: order of APA projection
%   - mu: step-size of the lineaar filter
%   - mQ: step-size of the nonlinearity
%   - delta: regularization parameter
%   - af: spline nonlinearity structure
%   - apa_af: structure of the created SAF filter
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


Mt = M + N;             % Total number of parametrs
b = zeros(M,1);         % MA filter taps
a = zeros(N,1);         % AR filter taps
w = zeros(Mt,1);        % Total AF taps (ARMA model)
beta  = zeros(M,N+1);   % Derivatives with respect MA coefficients
alpha = zeros(N,N+1);   % Derivatives with respect AR coefficients
dI = delta*eye(K);      % Regularizing diagonal matrix
X  = zeros(K,M);        % Covariance data matrix for APA
xd = zeros(K,1);        % Buffer of desired samples
xw = zeros(M,1);        % Buffer of filter status
sw = zeros(N,1);        % Buffer of IIR status
% -------------------------------------------------------------------------
apa_af = struct('M',M, 'N',N, 'K',K, 'mu',mu,'mQ',mQ, 'beta',beta, 'alpha',alpha, ...
    'dI',dI, 'b',b, 'a',a, 'w',w, 'X',X, 'xd',xd, 'xw',xw, 'sw',sw, 'af',af);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------