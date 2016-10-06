% WPOLY_demo.m
% 
% DEMO on nonlinear system identification (w0,q0) based on a Wiener
% Polynomial Adaptive Filter (WPOLY).
%
% This demo file implements the experimental test shown in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear Spline 
% Adaptive Filtering", Signal Processing, Vol. 93, No. 4, pp. 772-783, 
% April 2013.
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2012/05/15$
% $Revision: 1.1$  $Date: 2016/09/01$
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

clear all; 
close all;
disp('WPOLY_demo');
% -------------------------------------------------------------------------

%%  Parameters setting

% Input parameters --------------------------------------------------------
Lx = 50000;                  % Length of the input signal
nRun = 10;                   % Number of runs
out_noise_level_dB = 30;     % SNR
out_noise_level = 10^(-out_noise_level_dB/20);   % Noise level

x = zeros(Lx,1);             % x (Nx1) input signal array definition

% Colored signal generation -----------------------------------------------
a = 0.0;
b = sqrt(1-a^2);
% x = filter(b, [1 -a], randn( size(x)));   % H(z) = b/(1+a*z^-1)
disp('...... ');

% Adaptive filter definition ----------------------------------------------
M = 7;             % Length of the linear filter
K = 1;             % Affine Projection order; K=1 => NLMS
Pord = 3;          % Polynomial order
mu0 = 0.0005;      % Learning rate LMS/APA;
mQ0 = 0.0005;      % Learning rate for ctrl points;
if Lx < 30000,     % Batch for evaluating MSE
        B = 100;
else
        B = 4000;
end

% Spline activation function definition and initialization ----------------
afinit = 0;        % Init act. func. -1 0 ... (ONLY -1, bip.sig. or  0 =linear )
aftype = 3;        % Kind act. f. -1 0 1 2 3 4 5  
Slope  = 1;        % Slope
DeltaX = 0.2;      % Delta X
x_range = 2;       % Range limit


%% Initialization

% Creating the nonlinearity -----------------------------------------------
af0 = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M, Pord);   % Target
af1 = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M, Pord);   % Model

% --- Model definition ----------------------------------------------------
%H1 = create_SPL_apa_adaptive_filter_1(M,K,mu0,mQ0,1e-2,af1);   % Model APA
H1 = create_SPL_lms_adaptive_filter_1(M,mu0,mQ0,1e-2,af1);     % Model LMS

% --- Target Definition ---------------------------------------------------
%TH1 = create_SPL_apa_adaptive_filter_1(M,K,mu0,mQ0,1e-2,af0);   % Target Model APA
TH1 = create_SPL_lms_adaptive_filter_1(M,mu0,mQ0,1e-2,af0);     % Target Model LMS 

% --- Target Polynomial ---------------------------------------------------
Q0 = [ 1  -0.3  0.2 ].';   %3-th order
%Q0 = [ 1.4  0.7  1.2  0.9  1.2 ].';   %5-th order
%Q0 = [ 0.6984  0.1045  -0.0264  -0.0168  -0.0011  0.0008 ].';  %6-th order
%Q0 = [ 0.0008  -0.0011  -0.0168  -0.0264   0.1045   0.6984  -0.1968 ].';  %7-th order
%Q0 = [ -0.1968  0.6984  0.1045  -0.0264  -0.0168  -0.0011  0.0008 ].';  %7-th order
%Q0 = [ 1.2861  0.4493  -1.5836  -1.2253   1.4548  1.1157  -0.6093  -0.4541  0.1188  0.0871  -0.0089  -0.0065].';  % 12-th order
%Q0 = [1.7998  1.0069  -6.6350  -3.7890  15.9373  5.4953  -19.0208  -3.9872 12.1806  1.5391  -4.2507  -0.3003  0.7606  0.0232  -0.0546].'; % 15-th order
TH1.af.Q = Q0;
QL = length(Q0);
 
% --- Target linear part --------------------------------------------------
TH1.w  = [1  -0.4  0.25  -0.15  0.1  -0.05  0.001]';  % MA system to be identified

% Initialize --------------------------------------------------------------
N = Lx + M + K;
for i = Lx+1:N
    x(i)=0;
end

dn = zeros(N,1);        % Noise output array 
d  = zeros(N,1);        % Desired signal array
y  = zeros(N,1);        % Output array
e  = zeros(N,1);        % Error array 
em = zeros(N,1);        % Mean square error 
wm   = zeros(M,1);      % Mean value of w 
varW = zeros(M,1);      % Variance value of w 
qm   = zeros( QL, 1 );  % Mean value Spline coeff
varQ = zeros( QL, 1 );  % Variance value Spline coeff 


%% Main loop --------------------------------------------------------------
disp('WPOLY algorithm start ... ');
t = clock;

for n = 0 : nRun-1
    fprintf('Test nr. %d/%d\n', n+1, nRun);
    x = filter( b, [1 -a], randn( size(x)));  % H(z) = b/(1+a*z^-1)
    dn = out_noise_level*randn( size(x) );    % Noisy desired signal
       
    % WPOLY I.C. ----------------------------------------------------------
    H1.w (:) = 0;
    H1.w (1) = 1;    % Set filter i.c. 
        
    % Nonlinearity I.C. ---------------------------------------------------
    H1.af = create_activation_function(afinit, aftype,  DeltaX, x_range, Slope, M, Pord);
        
    % WPOLY Evaluation ----------------------------------------------------
    for k = 1 : Lx
        % Computing the desired output ------------------------------------
        [TH1, d(k), snk] = FW_WPOLY_F(TH1, x(k));       % Polynomial Wiener model     
        
        % Updating WPOLY --------------------------------------------------
        %[H1, y(k), e(k)] = AF_APA_WPOLY_F(H1, x(k), d(k) + dn(k) );  % WPOLY APA 
        [H1, y(k), e(k)] = AF_LMS_WPOLY_F(H1, x(k), d(k) + dn(k) );   % WPOLY LMS
    end

    em  = em  + (e.^2);   % Squared error

    % SAF run-time mean and variance estimation ---------------------------
    wm   = (1/(n+1))*H1.w + (n/(n+1))*wm;
    varW = varW + (n/(n+1))*((H1.w - wm).^2);
    
    qm   = (1/(n+1))*H1.af.Q + (n/(n+1))*qm;
    varQ = varQ + (n/(n+1))*((H1.af.Q - qm).^2); 
    
end

em = em/nRun;      % MSE
H1.af.Q = qm;

%--------------------------------------------------------------------------
% Average MSE evaluation
mse = mean(em(end-B-M-1:end-M-1));  % Average MSE
%--------------------------------------------------------------------------
fprintf('\n');


%% Results

% -------------------------------------------------------------------------
% Print table of means and variances
% -------------------------------------------------------------------------
fprintf('\n');
fprintf('Number of iterations = %d\n',nRun); 
fprintf('Learning rates:  muW = %5.3f   muQ = %5.3f\n', mu0, mQ0);
fprintf('a = %4.2f  b = %4.2f\n',a, b ); 
fprintf('Number of filter weights = %d\n', M); 
fprintf('Number of control points = %d\n', QL); 
fprintf('AF type = %d\n',aftype); 
fprintf('DeltaX  = %4.2f\n',DeltaX); 
fprintf('SNR_dB  = %4.2f dB\n',out_noise_level_dB);
fprintf('Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse,10*log10(mse));
fprintf('\n');
fprintf('Mean and Variance Tables -----------------------------------------\n');
for i=1:QL
    fprintf('i=%2d  q0 =%5.2f    qm =%9.6f   varQ = %10.3e \n', i,  H1.af.Q(i), qm(i), varQ(i) );
end
fprintf('\n');
fprintf('------------------------------------------------------------------\n');
for i=1:M
    fprintf('i=%d  w0 =%5.2f    wm =%9.6f   varW = %10.3e \n', i, H1.w(i), wm(i), varW(i) );
end
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
ylim([-yLIM yLIM]);   
xlim([-xLIM xLIM]);   
grid on;
KK = 500;
yy1 = zeros(1,KK);
yy2 = zeros(1,KK);
xa1 = zeros(1,KK);
dx = 2*xLIM/KK;
xx = -xLIM;   
for k = 1:KK
    yy1(k) = ActFunc(xx, TH1.af);   % Target  
    yy2(k) = ActFunc(xx, H1.af);    % WPOLY
    xa1(k) = xx;
    xx = xx + dx;
end  
title('Profile of polynomial nonlinearity after learning');
xlabel('Nonlinearity input {\itx}[{\itn}] ');
ylabel('Nonlinearity output {\its}[{\itn}]');
plot(xa1,yy1,'-.k','LineWidth',2);  % Target
plot(xa1,yy2,'r','LineWidth',2);    % WPOLY
legend('Target','WPLOY');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);

% MSE dB ------------------------------------------------------------------
figure2 = figure('PaperSize',[10 15]);
box('on');
hold on;
grid on;
ylim([-out_noise_level_dB-5 0]);   
[bb,aa] = butter(3, 0.02 );
edb = 10*log10( em ); 
plot( filter(bb,aa,edb ),'Color',[1 0 0],'LineWidth',2);   % MSE
noiseLevel(1: length(edb)-1 ) = -out_noise_level_dB;
plot( noiseLevel,'--','Color',[0 0 1],'LineWidth',2 );     % Noise level
title('Wiener Polynomial convergence test');
xlabel('Samples');
ylabel('MSE [dB]   10log({\itJ}({\bfw},{\bfQ}))');
legend('MSE', 'NoiseLevel' );
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);

% Filter coefficients -----------------------------------------------------
figure3 = figure('PaperSize',[20.98 29.68]);
hold on;
grid on;
title('Comparison of model and adapted linear filter');
plot(TH1.w,'--r','LineWidth',2.5);   % Target
plot(H1.w,'-k','LineWidth',2.5);     % WPLOY
xlabel('samples {\itn}');
ylabel('Linear combiner coefficients');
legend('Target','WPOLY');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);
% -------------------------------------------------------------------------

fprintf('END WPOLY_demo ---------------------------------------------------\n');

% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% University of Roma 'La Sapienza' 
% -------------------------------------------------------------------------