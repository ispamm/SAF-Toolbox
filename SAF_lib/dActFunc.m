function dx = dActFunc(y, af) 

% This function evaluate the derivative of the output of a spline 
% nonlinearity with respect its input.
%
% Details can be found in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear Spline 
% Adaptive Filtering", Signal Processing, Vol. 93, No. 4, pp. 772-783, 
% April 2013.
%
% USAGE:
%   dx = dActFunc(y, af)
%
% Input Arguments:
%   y       : input of the nonlinearity
%   af      : structure of the nonlinearity
%
% Output Argument:
%   dx      : derivative of the output of the nonlinearity
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

Sl = af.Slope;  % Slope 
G  = af.Gain;   % Gain 

% Nonlinearity selection
switch (af.aftype)
    
  case -1  % Signed Sigmoid 
      dx = ( Sl/(2*G)*(G^2-y^2) );  
      
  case  0  % Linear  
      dx = af.Slope ;   % Default value                      
      
  case  1  % Unsigned Sigmoid
      dx = ( (Sl/G)*y*(G-y) );          

  case  2   % Gaussian 
      dx = (-2*G/Sl)*af.s*exp( -af.s^2 / Sl );
  
  case  3   % Polynomial
      dx = af.Q(1);
      for k=2:af.Pord,
          dx = dx + k*af.Q(k)*y.^(k-1);
      end
    
  case 20   % Quadratic spline
      [~,af] = ActFunc(y,af);
      uIndex  = af.uIndex;    % Local index i
      u  = af.uSpline;        % Local abscissa u
      g  = [0 1 2*u]*af.C;    % Dot u vector         
      dx = g*af.Q( uIndex : uIndex+2 )/af.DeltaX;   % Eq. (6)
           
  otherwise   % Cubic spline  
      [~,af] = ActFunc(y,af);
      uIndex  = af.uIndex;       % Local iindex i
      u  = af.uSpline;           % Local abscissa u
      g  = [3*u^2 2*u 1 0]*af.C; % Dot u vector           
      dx = g*af.Q( uIndex : uIndex+3 )/af.DeltaX;   % Eq. (6)

end

% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------