% HSAF_demo.m
% 
% DEMO on nonlinear system identification (w0,q0) based on a Hammerstein 
% Spline Adaptive Filter (HSAF).
%
% This demo file implements the experimental test shown in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Hammerstein 
% Uniform Cubic Spline Adaptive Filters: Learning and Convergence 
% Properties", Signal Processing, Vol. 100, pp. 112-123, July 2014.
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2013/03/20$
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

clear all; 
close all;
disp('HSAF_demo');
% -------------------------------------------------------------------------

%%  Parameters setting

% Input parameters --------------------------------------------------------
Lx = 30000;               % Length of input signal
nRun = 10;                % Number of runs
out_noise_level_dB = 60;  % SNR
out_noise_level = 10^(-out_noise_level_dB/20);  % Noise level

x  = zeros(Lx,1);         % x (Nx1) input signal array definition

% Colored signal generation -----------------------------------------------
a = 0.1;
b = sqrt(1-a^2);
%x = filter( b, [1 -a], randn( size(x))) ;  % H(z) = b/(1+a*z^-1)
disp('...... ');

% Adaptive filter definition ----------------------------------------------
M =  7;           % Length of linear filter
mu0 = 0.1;        % Learning rate for linear filter
mQ0 = 0.1;        % Learning rate for control points
if Lx < 30000,    % Batch for evaluating MSE
        B = 100;
else
        B = 4000;
end

% Spline activation function definition and initialization ----------------
afinit = 0;    % Init act. func. -1 0 ... (ONLY -1, bip.sig. or  0 =linear ) 
aftype = 4;    % Kind act. f. -1 0 1 2 4 5; (4 = CR-spline, 5 = B-spline )
Slope  = 1;    % Slope
DeltaX = 0.2;  % Delta X
x_range = 2;   % Range limit

% Creating the nonlinearity -----------------------------------------------
af0 = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M);  % Model
af1 = create_activation_function(afinit, aftype, DeltaX, x_range, Slope, M);  % SAF


%% Initialization

% --- Target Definition ---------------------------------------------------
TH1 = create_SPL_lms_adaptive_filter_1(M,mu0,mQ0,1e-2,af0); % Target h Model LMS 

% TARGET: Nonlinear memoryless function implemented by Spline interpolated LUT
Q0 = [     -2.20 
           -2.00 
           -1.80 
           -1.60 
           -1.40 
           -1.20 
           -1.00 
           -0.80 
           -0.91 
           -0.40 
           -0.20 
            0.05 
            0.0  
           -0.40 
            0.58 
            1.00
            1.00 
            1.20 
            1.40 
            1.60 
            1.80 
            2.00
            2.20 
];
TH1.af.Q = Q0;
QL = length(Q0);   % Number of control points

% Linear filter -----------------------------------------------------------
TH1.w  = [0.6  -0.4  0.25  -0.15  0.1  -0.05  0.001]';  % MA system to be identified

% --- SAF definition ------------------------------------------------------
H1 = create_SPL_lms_adaptive_filter_1(M,mu0,mQ0,1e-2,af1); % HSAF LMS

% Initialize --------------------------------------------------------------
N = Lx + M + 1;   % Total samples
for i = Lx+1:N 
    x(i)=0;
end

dn = zeros(N,1);         % Noise desired output array
d  = zeros(N,1);         % Desired signal array
y  = zeros(N,1);         % Output array
e  = zeros(Lx,1);        % Error array 
em = zeros(Lx,1);        % Mean square error 
wm   = zeros(M,1);       % Mean value of w 
varW = zeros(M,1);       % Variance value of w 
qm   = zeros( QL, 1 );   % Mean value Spline coeff
varQ = zeros( QL, 1 );   % Variance value Spline coeff 


%% Main loop --------------------------------------------------------------
disp('Algorithm start ... ');
t = clock;

for n = 0 : nRun-1
    fprintf('Test nr. %d/%d\n', n+1, nRun); 
    x = filter( b, [1 -a], randn( size(x))) ;   % H(z) = b/(1+a*z^-1)
    dn = out_noise_level * randn( size(x) ) ;   % Noise
       
    % SAF I.C. ------------------------------------------------------------
    H1.w (:) = 0;
    H1.w (1) = 0.1;
    
    % Set Activation Func I.C. --------------------------------------------
    H1.af.Q = af0.Q;
    
    % HSAF Evaluation -----------------------------------------------------
    for k = 1 : Lx
        % Computing the desired output ------------------------------------
        [TH1, d(k), snk] = FW_HSPL_F(TH1, x(k));    % Hammerstein model
        
        % Updating HSAF ---------------------------------------------------
        [H1, y(k), e(k)] = AF_LMS_HSPL_F(H1, x(k), d(k) + dn(k) );  % SAF LMS ( Eqs. (10) and (11) ) 
    end

    em = em + (e.^2);   % Squared error
    
    % SAF run-time mean and variance estimation ---------------------------
    wm   = (1/(n+1))*H1.w + (n/(n+1))*wm;
    varW = varW  + (n/(n+1))*((TH1.w - wm).^2);
    
    qm   = (1/(n+1))*H1.af.Q + (n/(n+1))*qm;
    varQ = varQ     + (n/(n+1))*((TH1.af.Q - qm).^2);  
    
end
em  = em/nRun;  % MSE
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
    fprintf('i=%2d  q0 =%5.2f    qm =%9.6f   varQ = %10.3e \n', i,  TH1.af.Q(i), qm(i), varQ(i) );
end
fprintf('\n');
fprintf('------------------------------------------------------------------\n');
for i=1:M
    fprintf('i=%d  w0 =%5.2f    wm =%9.6f   varW = %10.3e \n', i, TH1.w(i), wm(i), varW(i) );
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
    yy1(k) = ActFunc(xx, TH1.af);   % Model  
    yy2(k) = ActFunc(xx, H1.af);    % Adapted
    xa1(k) = xx;
    xx = xx + dx;
end  
xlabel('Input {\itx}[{\itn}] ','FontSize', 12, 'FontWeight', 'demi');
ylabel('SAF output {\ity}[{\itn}]','FontSize', 12, 'FontWeight', 'demi');
title('Profile of model and adapted nonlinearity \varphi(x[n])','FontSize', 12, 'FontWeight', 'demi');
plot(xa1,yy1,'-.r','LineWidth',2);  
plot(xa1,yy2,'k','LineWidth',2);
legend('Target','Adapted','Location','SouthEast');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');

% Filter coefficients -----------------------------------------------------
figure2 = figure('PaperSize',[20.98 29.68]);
hold on;
plot(TH1.w,'LineWidth',2);
plot(H1.w,'m','LineWidth',2);
xlabel('time {\itn}','FontSize', 12, 'FontWeight', 'demi');
ylabel('Linear combiner coefficients','FontSize', 12, 'FontWeight', 'demi');
legend('Model','Adapted');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');

% MSE dB ------------------------------------------------------------------
figure3 = figure('PaperSize',[10 15]);
box('on');hold on;hold('all');
ylim([-out_noise_level_dB-5 10]);
grid on;
edb = 10*log10( em );
[bb,aa] = butter(2, 0.02 );
plot( filter(bb,aa,edb ),'Color',[1 0 0],'LineWidth',2);
noiseLevel(1: length(edb)-1 ) = -out_noise_level_dB;
plot( noiseLevel,'--','Color',[0 0 1],'LineWidth',2 );
title('Hammerstein SAF convergence test','FontSize', 12, 'FontWeight', 'demi');
xlabel('Samples','FontSize', 12, 'FontWeight', 'demi');
ylabel('MSE [dB]   10log({\itJ}({\bfw},{\bfQ}))','FontSize', 12, 'FontWeight', 'demi');
legend('MSE', 'NoiseLevel' );
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperSize', [20.98 29.68]);
% -------------------------------------------------------------------------

fprintf('END HSAF_demo ----------------------------------------------------\n');

% -------------------------------------------------------------------------
% © 2013
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% University of Roma 'La Sapienza' 
% -------------------------------------------------------------------------
