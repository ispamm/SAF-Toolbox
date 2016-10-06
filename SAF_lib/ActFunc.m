function [x, af] = ActFunc(s, af) 

% This function evaluate the output of a spline nonlinearity.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear Spline 
% Adaptive Filtering", Signal Processing, Vol. 93, No. 4, pp. 772-783, 
% April 2013.
%
% USAGE:
%   [x, af] = ActFunc(s, af)
%
% Input Arguments:
%   s       : input of the nonlinearity
%   af      : structure of the nonlinearity
%
% Output Arguments:
%   x       : output of the nonlinearity
%   af      : structure of the nonlinearity
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


af.s = s;    % Input of the nonlinearity

% Nonlinearity selection
switch (af.aftype)
    
  case -1  % Signed Sigmoidal
      x = (2*af.Gain/(1+exp(-s*af.Slope))-af.Gain);
  
  case  0  % Linear
      x = s*af.Slope; % Defaut value                   
      
  case  1  % Unsigned Sigmoidal
      x = af.Gain/(1+exp(-s*af.Slope));
      
  case  2  % Gaussian
      x = af.Gain*exp( -s^2/af.Slope );
      
  case  3  % Polynomial
      x = 0;
      for j=1:af.Pord,
          x = x + af.Q(j)*s.^j;  % Sum of monomials
          af.g(j) = s.^j;
      end
             
  case 20 % Quadratic spline
      np = af.lut_len;             % Number of control points
      Su = s/af.DeltaX + (np-1)/2; % First part of Eq. (7.b)
      uIndex = floor(Su);          % Span index i in Eq. (7.b)
      u = Su - uIndex;             % Local abscissa u in Eq. (7.a)
      if uIndex<1                  % The index must start from 1
          uIndex = 1;
      end
      if uIndex>(np-2),
          uIndex = np - 2;         % The index cannot exceed np - 2
      end
      af.g = [1 u u^2]*af.C;            % First part of Eq. (5): u^T C
      x = af.g*af.Q(uIndex : uIndex+2); % Eq. (5): u^T C q_i
      af.uIndex  = uIndex;              % For derivative computation
      af.uSpline = u;                   % For derivative computation   
      
  otherwise  % Cubic spline
      np = af.lut_len;             % Number of control points
      Su = s/af.DeltaX + (np-1)/2; % First part of Eq. (7.b)
      uIndex = floor(Su);          % Span index i in Eq. (7.b) 
      u = Su - uIndex;             % Local abscissa u in Eq. (7.a) 
      if uIndex<1                  % The index must start from 1
          uIndex = 1;
      end
      if uIndex>(np-3),            % The index cannot exceed np - 3
          uIndex = np - 3;
      end        
      af.g = [u^3 u^2 u 1]*af.C;        % First part of Eq. (5): u^T C
      x = af.g*af.Q(uIndex : uIndex+3); % Eq. (5): u^T C q_i
      af.uIndex  = uIndex;              % For derivative computation
      af.uSpline = u;                   % For derivative computation
   
end

% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------