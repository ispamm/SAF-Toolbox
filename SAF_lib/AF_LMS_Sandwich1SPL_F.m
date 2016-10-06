function [F, y, e] = AF_LMS_Sandwich1SPL_F(F, x, d)

% This function implements the adaptation of a Sandwich 1 SAF (S1SAF) structure 
% by using the LMS adaptive algorithm.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Novel Cascade 
% Spline Architectures for the Identification of Nonlinear Systems", 
% IEEE Transactions on Circuits and Systems---I: Regular Papers, Vol. 62,
% No. 7, pp. 1825-1835, July 2015.
%
% USAGE:
%  [F, y, e] =  AF_LMS_Sandwich1SPL_F(F, x, d)
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

M = F.M;         % Length of the linear filter

F.xw(2 : M) = F.xw(1 : M-1);           % Shift the input delay-line
[F.xw(1), F.af1] = ActFunc(x, F.af1);  % Load a new input into the delay-line

for j=2:M,          % Construct the U matrix
    if F.af1.uIndex == F.af1.indexes(j)
        F.af1.gM(j,:) = F.af1.gM(j-1,:);
    else
        F.af1.gM(j,:) = zeros(1,4);
    end
end
F.af1.gM(1,:) = F.af1.g;
F.af1.indexes(2:M) = F.af1.indexes(1:M-1);
F.af1.indexes(1) = F.af1.uIndex;

r = F.xw.'*F.w;                % Linear filter output array 
[y,F.af2] = ActFunc(r,F.af2);  % Output array
e = d - y;                     % Error array
ee = e*dActFunc(r,F.af2);      % Error multiplied the nonlinearity derivative

% LMS weights and control points update -----------------------------------
F.w = F.w + F.mu*conj(ee)*F.xw;   % cpx LMS in Eq. (30)

if F.af1.aftype > 1  % Kind act. f. -1 0 1 2 4 5 *
    e_av = F.mQ1*ee;     
    ii = F.af1.uIndex : F.af1.uIndex + F.af1.P;  % P = Spline order

    F.af1.Q(ii) = F.af1.Q(ii) + e_av*(F.w'*F.af1.gM).';  % LMS in Eq. (31)
end 

if F.af2.aftype > 1  % Kind act. f. -1 0 1 2 4 5 *
    e_av = F.mQ2*e;     
    ii = F.af2.uIndex : F.af2.uIndex + F.af2.P;  % P = Spline order

    F.af2.Q(ii) = F.af2.Q(ii) + e_av*F.af2.g.';   % LMS in Eq. (32)
end 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2014
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------