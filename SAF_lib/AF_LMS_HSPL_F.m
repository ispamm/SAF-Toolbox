function  [F, y, e] =  AF_LMS_HSPL_F(F, x, d)

% This function implements the adaptation of a Hammerstein SAF (HSAF) structure 
% by using the LMS adaptive algorithm.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Hammerstein 
% Uniform Cubic Spline Adaptive Filters: Learning and Convergence 
% Properties", Signal Processing, Vol. 100, pp. 112-123, July 2014.
%
% USAGE:
%  [F, y, e] =  AF_LMS_HSPL_F(F, x, d)
%
% Input Arguments:
%   F       : HSAF structurs
%   x       : array of input signal
%   d       : array of desired signal
%  
% Output Arguments:
%   F       : HSAF structure
%   y       : array of output signal   
%   e       : array of error signal    
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2013/05/17$
% $Revision: 1.1$  $Date: 2016/07/31$
% -------------------------------------------------------------------------
% License to use and modify this code is granted freely without warranty to
% all, as long as the original authors are referenced and attributed as such.
% The original authors maintain the right to be solely associated with this
% work.
% -------------------------------------------------------------------------
% © 2013
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------


M = F.M;                             % Length of the linear filter
F.xw(2 : M) = F.xw(1 : M-1);         % Shift the input delay-line
[F.xw(1), F.af] = ActFunc(x, F.af);  % Load a new input into the delay-line  
for j=2:M,                           % Constructing the matrix U
    if F.af.uIndex == F.af.indexes(j)
        F.af.gM(j,:) = F.af.gM(j-1,:);
    else
        F.af.gM(j,:) = zeros(1,4);
    end
end
F.af.gM(1,:) = F.af.g;               % Load a new vector in matrix U
F.af.indexes(2:M) = F.af.indexes(1:M-1);
F.af.indexes(1) = F.af.uIndex;

y = F.xw.'*F.w;                      % Filter output 
e = d - y;                           % Error estimation


% LMS weights and control points update -----------------------------------
F.w = F.w + F.mu*F.xw.*e;            % LMS in Eq. (10)

if F.af.aftype > 1  % Kind act. f. -1 0 1 2 4 5 *
    e_av = F.mQ*e;     
    ii = F.af.uIndex : F.af.uIndex + F.af.P;         % P = Spline order
    F.af.Q(ii) = F.af.Q(ii) + e_av*(F.w'*F.af.gM).'; % LMS in Eq. (11)
end 
%*************************************************************************

% -------------------------------------------------------------------------
% © 2013
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------