function AF = create_activation_function(afinit, aftype, DeltaX, Gain, Slope, M, Pord)

% Spline nonlinear function creation for SAF.

% Details can be found in:
% [1] M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear 
%     Spline Adaptive Filtering", Signal Processing, Vol. 93, No. 4, 
%     pp. 772-783, April 2013.
%
% [2] M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Hammerstein 
%     Uniform Cubic Spline Adaptive Filters: Learning and Convergence 
%     Properties", Signal Processing, Vol. 100, pp. 112-123, July 2014.
%
%
% USAGE:
%    AF = create_activation_function(afinit, aftype, DeltaX, Gain, Slope, M, Pord)
% 
% Paramaters:
%    afinit  : Shape of init activation function (i.c. if flexible) 
%                    -1 signed sigmoid 
%                     0 linear function
%                     1 unsigned sigmoid
%                     2 Gaussian
%                     3 random
%    aftype  : Kind of activation function -1 0 1 2 4 5  
%                    -1 fixed signed sigmoid 
%                     0 fixed linear function
%                     1 fixed unsigned sigmoid
%                     2 fixed Gaussian
%                     3 fexible polynomial
%                     4 fexible Catmull-Rom spline
%                     5 fexible B-spline
%   DeltaX   : X-axes sampling interval for flexible a.f.
%   Slope    : Slope of a.f.
%   Gain     : Amplitude of a.f.
%   M        : Length of the linear filter
%   Pord     : Order of polynomial (for polynomial nonlinearity)
%   AF       : structure containing the nonlinear block
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



% Spline activation function definition and initialization ----------------
if nargin==0, help create_activation_function; return; end
if nargin < 7
    Pord = 3;
    if nargin<6
        M=1;
        if nargin<5 
            Slope=2.16;  
            if nargin<4 
                Gain=1.1;
                if nargin<3
                    DeltaX=0.4;
                    if nargin<2
                        aftype=-1;
                        if nargin<1
                            afinit=-1;
                        end
                    end
                end
            end 
        end
    end
end

% -------------------------------------------------------------------------
% Check fo nonlinearity type
if afinit == -1, 
    Slope  = Slope*2*(1/Gain);
end

if afinit == 1,
    Slope  = Slope*4*(1/Gain);
end

if ((afinit<-1) || (afinit>3) )
    fprintf('Activation function error type\n'); 
    aftype = -1;
end

% LUT parameters ----------------------------------------------------------
Table_Length = TabFuncLen(DeltaX, Gain ,Slope, afinit );
lut_len      = Table_Length; % Length of LUT
s = 0.0;     % Linear combiner output
x = 0.0;     % Nonlinearity output
uIndex = 1;  % Span index
uSpline = 0; % Local abscissa

% Polynomial nonlinearity -------------------------------------------------
if (aftype == 3)
    Q = zeros(Pord,1);     % Polynomial nonlinearity
    Q(1) = 1;
    g = zeros(1,Pord);     % Dot u vector
    gM = zeros(Pord,M);    % U matrix (for Hammerstein filter) in [2]
else
    Q = zeros(lut_len,1);  % Look-Up Table nonlinearity
end

% If NOT spline -----------------------------------------------------------
if (aftype < 4) 
    C = 0;
end; 

% Catmul-Rom spline nonllinearity -----------------------------------------
if (aftype == 4)
C = 0.5*[-1  3 -3  1; ...
          2 -5  4 -1; ...  
         -1  0  1  0; ...          % Page 775, top of column 1
          0  2  0  0];
end

% B spline nonllinearity --------------------------------------------------
if (aftype == 5)
    C = (1/6)*[-1  3 -3  1; ... 
                3 -6  3  0; ...     % Page 774, bottom of column 2
               -3  0  3  0; ...  
                1  4  1  0];
end

% Bernstein polynomial nonlinearity ---------------------------------------
if (aftype == 6)
    C =       [-1  3 -3  1; ... 
                3 -6  3  0; ...   % Bernstein polynomials
               -3  3  0  0; ...  
                1  0  0  0];
end

% Parametric spline nonlinearity ------------------------------------------
if (aftype == 7)
     tau = 0.5;   % tau = 0.5 ==> CR-spline
     C =       [-tau   2-tau  tau-2   tau; ... 
                 2*tau tau-3  3-2*tau -tau; ...   % Parametric spline
                -tau   0      tau      0 ; ...  
                 0     1      0        0 ];
end

% Hermite spline nonlinearity ---------------------------------------------
if (aftype == 8)
    C =       [ 2 -2  1  1; ... 
               -3  3 -2 -1; ...   % Hermite polynomials
                0  0  1  0; ...  
                1  0  0  0];
end

% Bezier spline nonllinearity ---------------------------------------------
if (aftype == 9) 
    C = (1/6)*[ 1  3 -3  1; ... 
                3 -6  3  0; ...   % Bezier polynomials
               -3  3  0  0; ...  
                1  0  0  0];
end

% Quadratic B-spline nonlinearity -----------------------------------------
if (aftype == 20) 
     C =       [ 1  1  0; ... 
                -2  2  0; ...
                 1 -2  1];
end

% Spline order ------------------------------------------------------------
P = length(C) - 1;

% Dummy arrays for learning algorithms. See Eqs. (5) and (6) --------------
if (aftype > 3)
    g = zeros(1,P+1);   % The dot u vector
    gM = zeros(M,P+1);  % The U matrix (for Hammerstein SAF) in [2]
end
indexes = ones(M,1);

% STRUCTURE DEFINITION ----------------------------------------------------
AF = struct('aftype', aftype, 'afinit',afinit, 's',s, 'x',x, 'Q',Q, ...
            'lut_len', lut_len, 'uIndex',uIndex, 'uSpline',uSpline, ...
            'Slope',Slope, 'Gain',Gain, 'DeltaX',DeltaX, ...
            'P', P, 'Pord', Pord, 'C',C, 'g',g, 'gM',gM, 'indexes',indexes);
% -------------------------------------------------------------------------

% For spline interpolation ------------------------------------------------
if (aftype>1 && aftype ~= 3)
    LutSlope = (Table_Length - 1)/2.0;  % New slope 
    X = -LutSlope*DeltaX; 
    for j=1 : Table_Length   % Table_Length   
        AF.Q(j) = FUNC(X, Gain, Slope, afinit);
        X = X + DeltaX;
    end 
end
% -------------------------------------------------------------------------


%=============== Function TabFuncLen() ====================================
function Table_Lenght = TabFuncLen(DX, G, S, ty)
    ii = 0;
    X  = 0;
    if (ty~=2 && ty~=3)   
        F = 0.0;
        crtGain = G - 0.005*G; % Max atc func value -----------------------
        while (F<crtGain)
            F  = FUNC(X,G,S,ty);
            X  = X + DX;
            ii = ii + 1;      
        end
    elseif (ty==3)
        ii = 11;
    else
        crtGain = 0.005*G;  % Gaussian
        F = G;
        while (F>crtGain )
             ii = ii + 1;          
             F  = FUNC(X,G,S,ty);
             X  = X + DX;
        end
    end
    Table_Lenght = ii*2 + 1;   % always odd
% End function TabFuncLen() ===============================================

% =================== Function FUNC() =====================================
function value = FUNC(X, G, S, ty)
% Nonlinear function implementation
    switch  ty
        case -1  
            value = 2*G/(1 + exp(-X*S) ) - G;  % Signed Sigmoid
        case  0  
            value = X*S;                % Linear
        case  1  
            value = G/(1+exp(-X*S));    % Unsigned Sigmoidal
        case  2  
            value = G*exp(-(X*X)/5.0);  % Gaussian
        case 3
            value = 2*G*(rand - 0.5);   % Random
        otherwise
            value = 0;
    end
% End function FUNC() =====================================================

% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------