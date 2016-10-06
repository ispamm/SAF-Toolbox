% HSAF_Polynomial_Compare_demo.m
% 
% DEMO of comparison results on nonlinear system identification based on a 
% Hammerstein Spline Filter (HSAF) with a III-order Volterra Adaptive Filter
% and Polynomial Adaptive Filters.
%
% This demo file implements the experimental test shown in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Hammerstein 
% Uniform Cubic Spline Adaptive Filters: Learning and Convergence 
% Properties", Signal Processing, Vol. 100, pp. 112-123, July 2014.
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2013/03/20$
% $Revision: 1.1$  $Date: 2016/09/03$
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

clear all; 
close all;
disp('HSAF_Polynomial_Compare_demo');
% -------------------------------------------------------------------------

%%  Parameters setting

% Input parameters --------------------------------------------------------
Lx = 50000;                    % Length of the input signal
nRun = 10;                     % Number of runs
out_noise_level_dB = Inf;      % SNR
out_noise_level = 10^(-out_noise_level_dB/20);   % Noise level

x = zeros(Lx,1);               % x (Nx1) input signal array definition

% Colored signal generation -----------------------------------------------
a = 0.0;
b = sqrt(1-a^2);
% x = filter(b, [1 -a], randn( size(x)));   % H(z) = b/(1+a*z^-1)
disp('...... ');

% Adaptive filter definition ----------------------------------------------
M  = 6;             % Length of the linear part of HSAF
MV = 15;            % Length of the Volterra filter
K  =  1;            % Affine Projection order; K=1 => NLMS
Kv =  1;            % Volterra Affine Projection order; K=1 => NLMS
delta = 1e-2;       % APA regularization parameters
mu  = 0.01;         % Learning rate for the linear filter of HSAF
mQ  = 0.02;         % Learning rate for ctrl points;
mu2 = 0.001;        % Learning rate for the linear part of Polynomial filter
mQ2 = 0.001;        % Learning rate for the nonlinear part
mu3 = 0.00001;      % Learning rate for the linear filter of Mathews system
mQ3 = 0.00001;      % Learning rate for the nonlinear filter
muV = 0.01;         % Learning rate for the Volterra filter
if Lx < 30000,     % Batch for evaluating MSE
        B = 100;
else
        B = 4000;
end

% Spline activation function definition and initialization ----------------
afinit  = 0;        % Kind of init act. f. -1 0 1 2 3 */
aftype  = 5;        % Kind act. f. -1 0 1 2 3 4 5
Slope   = 1;        % Slope
DeltaX  = 0.2;      % Delta X
x_range = 2;        % Range limit
aftype2 = 3;        % Kind act. f. -1 0 1 2 3 4 5
Pord    = 3;        % Polynomial Order


%% Initialization

% Creating the nonlinearity -----------------------------------------------

% Targets -----------------------------------------------------------------
[af01]  = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M);         % HSAF
[af02]  = create_activation_function(afinit, aftype2, DeltaX, x_range, Slope, M, Pord);  % Polynomial
[af03]  = create_activation_function(afinit, aftype2, DeltaX, x_range, Slope, M, Pord);  % Mathews

% Models ------------------------------------------------------------------
[af1]   = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M);         % HSAF
[af2]   = create_activation_function(afinit, aftype2, DeltaX, x_range, Slope, M, Pord);  % Polynomial
[af3]   = create_activation_function(afinit, aftype2, DeltaX, x_range, Slope, M, Pord);  % Mathews

% --- Models definition ---------------------------------------------------
apa_af1 = create_SPL_apa_adaptive_filter_1(M, K, mu, mQ, delta, af1);   % HSAF APA
%apa_af1 = create_SPL_lms_adaptive_filter_1(M, mu, mQ, delta, af1);      % HSAF LMS
apa_af2 = create_SPL_apa_adaptive_filter_1(M, K, mu2, mQ2, delta, af2); % Polynomial APA   
%apa_af2 = create_SPL_lms_adaptive_filter_1(M, mu2, mQ2, delta, af2);    % Polynomial LMS
apa_af3 = create_SPL_lms_adaptive_filter_1(M, mu2, mQ2, delta, af3);    % Mathews LMS
%VF    = create_II_ord_Volterra_filter_1(MV, Kv, muV, delta);            % II Order Volterra
VF    = create_III_ord_Volterra_filter_1(MV, Kv, muV, delta);           % III Order Volterra


% Initialize --------------------------------------------------------------
N = Lx + M + K;
for i = Lx+1:N
    x(i)=0;
end

dn = zeros(N,1);        % Noise output array 
d  = zeros(N,1);        % Desired signal array
y1  = zeros(N,1);       % Output array 
y2  = zeros(N,1);       % Output array 
y3  = zeros(N,1);       % Output array 
y4  = zeros(N,1);       % Output array
e1  = zeros(Lx,1);      % Error array
e2  = zeros(Lx,1);      % Error array 
e3  = zeros(Lx,1);      % Error array 
e4  = zeros(Lx,1);      % Error array
em1 = zeros(Lx,1);      % Mean square error 
em2 = zeros(Lx,1);      % Mean square error
em3 = zeros(Lx,1);      % Mean square error
em4 = zeros(Lx,1);      % Mean square error


%% Main loop --------------------------------------------------------------
disp('HSAF_Polynomial_compare algorithm start ... ');
t = clock;

for n = 0:nRun-1
    fprintf('Test nr. %d/%d\n', n+1, nRun);
    
    % Nonlinear system signal generation ----------------------------------
    x = randn(Lx,1);           % x (Nx1) input signal array definition
    if a>0,
        x = filter( b, [1 -a], x);               % Colored input
    end
    z = 1/8 + fix((8*x - 4)/4);                  % Nonlinearity
    h = zeros(1,6);
    for k=0:5
        h(k+1) = 1 - k/4;
    end        
    d = 0.5*filter(h, 1, z);                     % Desired signal
    dn = out_noise_level*randn( size(x) );       % Noisy desired signal
    
    % Models I.C. ---------------------------------------------------------
    apa_af1.w(:) = 0;
    apa_af1.w(1) = 1;
    apa_af2.w(:) = 0;
    apa_af2.w(1) = 1;
    apa_af3.w(:) = 0;
    apa_af3.w(1) = 1;
    %VF    = create_II_ord_Volterra_filter_1(MV, Kv, muV, delta);   % II Order Volterra
    VF    = create_III_ord_Volterra_filter_1(MV, Kv, muV, delta);  % III Order Volterra
    
    % Nonlinearities I.C. -------------------------------------------------
    af1.Q  = af01.Q;
    af2.Q  = af02.Q;
    af3.Q  = af03.Q;
    
   % Models Evaluation ----------------------------------------------------
    for k = 1 : Lx  

        % Updating HSAF ---------------------------------------------------
        [apa_af1, y1(k), e1(k)] = AF_APA_HSPL_F(apa_af1, x(k), d(k) + dn(k) );          % HSAF
        %[apa_af1, y1(k), e1(k)] = AF_LMS_HSPL_F(apa_af1, x(k), d(k) + dn(k) );          % HSAF
        
        % Updating HPOLY --------------------------------------------------
        [apa_af2, y2(k), e2(k)] = AF_APA_HPOLY_F(apa_af2, x(k), d(k) + dn(k) );         % HPOLY
        %[apa_af2, y2(k), e2(k)] = AF_LMS_HPOLY_F(apa_af2, x(k), d(k) + dn(k) );         % HPOLY
               
        % Updating Mathews ------------------------------------------------
        [apa_af3, y3(k), e3(k)] = AF_LMS_MATHEWS_POLY_F(apa_af3, x(k), d(k) + dn(k) );  % MATHEWS 
                
        % Updating VAF ----------------------------------------------------
        %[VF, y4(k), e4(k)] =  AF_II_ord_VF_APA_F( VF, x(k), d(k) + dn(k) );             % II Order Volterra
        [VF, y4(k), e4(k)] =  AF_III_ord_VF_APA_F( VF, x(k), d(k) + dn(k) );            % III Order Volterra
    end
    
    em1  = em1  + (e1.^2);   % Squared error
    em2  = em2  + (e2.^2);   % Squared error
    em3  = em3  + (e3.^2);   % Squared error
    em4  = em4  + (e4.^2);   % Squared error
    
end

em1 = em1/nRun;      % MSE
em2 = em2/nRun;      % MSE
em3 = em3/nRun;      % MSE
em4 = em4/nRun;      % MSE

%--------------------------------------------------------------------------
% Average MSE evaluation
mse1 = mean(em1(end-B-M-1:end-M-1));  % Average MSE
mse2 = mean(em2(end-B-M-1:end-M-1));  % Average MSE
mse3 = mean(em3(end-B-M-1:end-M-1));  % Average MSE
mse4 = mean(em4(end-B-M-1:end-M-1));  % Average MSE
%--------------------------------------------------------------------------
fprintf('\n');


%% Results

% -------------------------------------------------------------------------
% Print table
% -------------------------------------------------------------------------

fprintf('\n');
fprintf('Mean Square Errors -----------------------------------------------\n');
fprintf('      HSAF Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse1,10*log10(mse1));
fprintf('Polynomial Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse2,10*log10(mse2));
fprintf('   Mathews Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse3,10*log10(mse3));
fprintf('       VAF Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse4,10*log10(mse4));
fprintf('------------------------------------------------------------------\n');


% -------------------------------------------------------------------------
% Plotting figures
% -------------------------------------------------------------------------

% Plot Spline functions ---------------------------------------------------
yLIM = 1.5;
xLIM = 3.0;
figure1 = figure('PaperSize',[20.98 29.68]);
box('on');
hold on;
hold('all');
grid on;
ylim([-yLIM yLIM]);   
xlim([-xLIM xLIM]);   
KK = 500;
yy1 = zeros(1,KK);
yy2 = zeros(1,KK);
yy3 = zeros(1,KK);
yy4 = zeros(1,KK);
yy5 = zeros(1,KK);
xa1 = zeros(1,KK);
dx = 2*xLIM/KK;
xx = -xLIM;   
for k = 1:KK
    yy1(k) = ActFunc(xx, apa_af1.af);  % Spline   
    yy2(k) = ActFunc(xx, apa_af2.af);  % Polynomial
    yy3(k) = ActFunc(xx, apa_af3.af);  % Mathews
    [VF,yy4(k)] = FW_III_ord_VF_F(VF, xx); % Volterra
    %yy5(k) = real(xx.^0.33333333);     % Target 33
    %yy5(k) = nthroot(xx,3);            % Target 33
    yy5(k) = 1/8 + fix((8*xx - 4)/4);  % Target 34
    %yy5(k) = xx - 0.3*xx.^2 + 0.2*xx.^3; % Target 103
    xa1(k) = xx;
    xx = xx + dx;
end  
title('Profile of nonlinearities after learning','FontSize', 12, 'FontWeight', 'demi');
xlabel('Nonlinearity input {\itx}[{\itn}] ','FontSize', 12, 'FontWeight', 'demi');
ylabel('Nonlinearity output {\its}[{\itn}]','FontSize', 12, 'FontWeight', 'demi');
plot(xa1,yy1,'k','LineWidth',2);     % HSAF
plot(xa1,yy2,'g','LineWidth',2);     % Polynomial
plot(xa1,yy3,'r','LineWidth',2);     % Mathews
plot(xa1,yy4,'c','LineWidth',2);     % Volterra
plot(xa1,yy5,'-.b','LineWidth',2);   % Target
legend('HSAF','Polynomial','Mathews','Volterra','Location','SouthEast');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);

% MSE dB ------------------------------------------------------------------
figure2 = figure('PaperSize',[20.98 29.68]);
box('on');
hold on;
grid on;
ylim([-out_noise_level_dB-5 0]);   
[bb,aa] = butter(2, 0.005);
edb1 = 10*log10( em1 );
edb2 = 10*log10( em2 );
edb3 = 10*log10( em3 );
edb4 = 10*log10( em4 );
plot(filter(bb,aa,edb1),'Color',[0 0 0],'LineWidth',2);       % HSAF
plot(filter(bb,aa,edb2),'Color',[0 1 0],'LineWidth',2);       % Poly 
plot(filter(bb,aa,edb3),'Color',[1 0 0],'LineWidth',2);       % Mathews
plot(filter(bb,aa,edb4),'Color',[0.2 0.9 0.8],'LineWidth',2); % Volterra
title('Comparisons of MSE','FontSize', 12, 'FontWeight', 'demi');
xlabel('Samples','FontSize', 12, 'FontWeight', 'demi');
ylabel('MSE [dB]   10log({\itJ}({\bfw},{\bfQ}))','FontSize', 12, 'FontWeight', 'demi');
legend('HSAF','Polynomial','Mathews','3-rd Order Volterra');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2013
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------