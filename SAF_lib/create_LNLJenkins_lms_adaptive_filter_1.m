function lms_af = create_LNLJenkins_lms_adaptive_filter_1(M1, M2, M3, mu1, mu2, mu3, delta, af)

% This function creates and initializes the structure implementing an
% adaptive LNL sandwich model proposed in:
% V. Hegde, C. Radhakrishnan, D. J. Krusienski, & W. K. Jenkins, (2002), 
% Series-Cascade Nonlinear Adaptive Filters, in 'Proceedings of the 45th 
% Midwest Symposium on Circuits and Systems (MWSCAS2002)', pp. 219--222.
%
% It has been used for comparisons in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Novel Cascade 
% Spline Architectures for the Identification of Nonlinear Systems", 
% IEEE Transactions on Circuits and Systems---I: Regular Papers, Vol. 62,
% No. 7, pp. 1825-1835, July 2015.
%
% USAGE:
% lms_af = create_LNLJenkins_lms_adaptive_filter_1(M1, M2, M3, mu1, mu2, mu3, delta, af))
%   - M1: length of the first part of LNL filter
%   - M2: length of the second part of LNL filter
%   - M3: length of the third part of LNL filter
%   - mu1: step-size of the first part of LNL filter
%   - mu2: step-size of the second part of LNL filter
%   - mu3: step-size of the third part of LNL filter
%   - delta: regularization parameter
%   - af: spline nonlinearity structure
%   - lms_af: structure of the created LNL filter
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


% LNL Jenkins adaptive filter definition and initialization ---------------
w1  = zeros(M1,1);       % Taps of the first part of LNL filter  
w2  = zeros(M2,1);       % Taps of the first second of LNL filter  
w3  = zeros(M3,1);       % Taps of the first third of LNL filter
w1(floor(M1/2)) = 1;     % Adaptive filter weights i.c.
w2(floor(M2/2)) = 1;     % Adaptive filter weights i.c.
dI  = delta;             % Regularizing diagonal matrix
xd  = zeros(M3,1);       % Buffer of the desired samples
xw1 = zeros(M1,1);       % Buffer of the filter 1 status
xw2 = zeros(M2,1);       % Buffer of the filter 2 status
xw3 = zeros(M3,1);       % Buffer of the filter 3 status
X   = zeros(M1,M3);      % Buffer matrix of the past inputs delay lines
Z   = zeros(M2,M3);      % Buffer matrix of the past nonlinear delay lines
% -------------------------------------------------------------------------
lms_af = struct('M1',M1, 'M2',M2, 'M3',M3, 'mu1',mu1, 'mu2',mu2, 'mu3',mu3, ...
    'dI',dI, 'w1',w1, 'w2',w2, 'w3',w3, 'xd',xd, 'xw1',xw1, 'xw2',xw2, ... 
    'xw3',xw3, 'X',X, 'Z',Z, 'af',af);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------