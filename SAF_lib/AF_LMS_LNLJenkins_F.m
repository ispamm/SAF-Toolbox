function [F, y, e] = AF_LMS_LNLJenkins_F(F, x, d)

% This function implements the adaptive LNL sandwich model proposed in:
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
%  [F, y, e] =  AF_LMS_LNLJenkins_F(F, x, d)
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


M1 = F.M1;                        % Length of the first filter
M2 = F.M2;                        % Length of the second filter
M3 = F.M3;                        % Length of the thirs filter

F.xw1(2 : M1) = F.xw1(1 : M1-1);  % Shift the filter 1 delay-line
F.xw1(1) = x;                     % Load a new input in the delay-line

F.xw2(2 : M2) = F.xw2(1 : M2-1);  % Shift the filter 2 delay-line
F.xw2(1) = F.xw1'*F.w1;           % Load a new sample in the delay-line
z = zeros(F.af.Pord,1);
for j=1:F.af.Pord,
    z(j) = F.xw2(1).^j;           % Evaluate the polynomial nonlinearity
end

F.xw3(2 : M3) = F.xw3(1 : M3-1);  % Shift the filter 3 delay-line
F.xw3(1) = z'*F.w2;               % Load a new sample in the delay-line

y = F.xw3'*F.w3;                  % Output evaluation
e = d - y;                        % Error evaluation

zp = zeros(F.af.Pord,1);
for j=2:F.af.Pord,
    zp(j) = j*F.xw2(1).^(j-1);    % Evaluate the polynomial nonlinearity
end
zp(1) = 1;
zq = zp'*F.w2;
F.X(:,2:M3) = F.X(:,1:M3-1);
F.X(:,1) = zq.*F.xw1;
F.Z(:,2:M3) = F.Z(:,1:M3-1);
F.Z(:,1) = z;

% LMS weights update ------------------------------------------------------
F.w1 = F.w1 + F.mu1*conj(e)*F.X*F.w3;   % Updatade first filter
F.w2 = F.w2 + F.mu2*conj(e)*F.Z*F.w3;   % Updatade second filter
F.w3 = F.w3 + F.mu3*conj(e)*F.xw3;      % Updatade third filter

% NLMS weights update -----------------------------------------------------
% F.w1 = F.w1 + F.mu1*conj(e)*F.X*F.w3/(F.dI + F.w3'*F.X'*F.X*F.w3);
% F.w2 = F.w2 + F.mu2*conj(e)*F.Z*F.w3/(F.dI + F.w3'*F.Z'*F.Z*F.w3);
% F.w3 = F.w3 + F.mu3*conj(e)*F.xw3/(F.dI + F.xw3'*F.xw3);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------