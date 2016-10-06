function [F, y, e] = AF_LMS_Sandwich2SPL_F(F, x, d)

% This function implements the adaptation of a Sandwich 2 SAF (S2SAF) structure 
% by using the LMS adaptive algorithm.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Novel Cascade 
% Spline Architectures for the Identification of Nonlinear Systems", 
% IEEE Transactions on Circuits and Systems---I: Regular Papers, Vol. 62,
% No. 7, pp. 1825-1835, July 2015.
%
% USAGE:
%  [F, y, e] =  AF_LMS_Sandwich2SPL_F(F, x, d)
%
% Input Arguments:
%   F       : WSAF structurs
%   x       : array of input signal
%   d       : array of desired signal
%  
% Output Arguments:
%   F       : WSAF structure
%   y       : array of output signal   
%   e       : array of error signal  
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


M1 = F.M1;                        % Length of the first linear filter
M2 = F.M2;                        % Length of the seconf linear filter

F.xw1(2 : M1) = F.xw1(1 : M1-1);  % Shift the input filter 1 delay-line
F.xw1(1) = x;                     % Load a new input into the delay-line
F.xs(2 : M2)  = F.xs(1 : M2-1);   % Shift the nonlinearity delay-line
F.xs(1) = F.xw1'*F.w1;            % Load a new input into the delay-line  
F.xw2(2 : M2) = F.xw2(1 : M2-1);  % Shift the filter 2 delay line
[F.xw2(1),F.af] = ActFunc(F.xs(1),F.af);  % Load a new input into the delay-line  
F.X(:,2:M2) = F.X(:,1:M2-1);      % Shift the data matrix delay-line
F.X(:,1) = F.xw1;                 % Load a new input into the matrix delay-line  
y = F.xw2'*F.w2;                  % Output array
e = d - y;                        % Error array
ee = e*dActFunc(F.xs(1),F.af);    % Error multiplied the nonlinearity derivative

for j=2:M2,   % Update the U matrix
    if (F.af.uIndex >= F.af.indexes(j)-F.af.P) && (F.af.uIndex <= F.af.indexes(j)+F.af.P)
        F.af.gM(j,:) = F.af.gM(j-1,:);
    else
        F.af.gM(j,:) = zeros(1,4);
    end
end
F.af.gM(1,:) = F.af.g;
F.af.indexes(2:M2) = F.af.indexes(1:M2-1);
F.af.indexes(1) = F.af.uIndex;

% LMS weights and LUT control points update -------------------------------
F.w1 = F.w1 + F.mu1*conj(ee)*F.X*F.w2;  % cpx LMS in Eq. (32)
F.w2 = F.w2 + F.mu2*conj(e)*F.xw2;      % cpx LMS in Eq. (30)

if F.af.aftype > 1  % Kind act. f. -1 0 1 2 4 5 *
    e_av = F.mQ*e;     
    ii = F.af.uIndex : F.af.uIndex + F.af.P;  % P = Spline order

    F.af.Q(ii) = F.af.Q(ii) + e_av*F.af.gM.'*F.w2;   % LMS in Eq. (31)
end 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------