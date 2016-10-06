% WSAF_DX_compare_demo.m
% 
% DEMO on nonlinear system identification (w0,q0) based on a Wiener 
% Spline Adaptive Filter (WSAF) using different values of Delta X.
%
% This demo file implements the experimental test shown in:
% M. Scarpiniti, D. Comminiello, R. Parisi and A. Uncini, "Nonlinear Spline 
% Adaptive Filtering", Signal Processing, Vol. 93, No. 4, pp. 772-783, 
% April 2013.
%
% Info: michele.scarpiniti@uniroma1.it
% -------------------------------------------------------------------------
% $Revision: 1.0$  $Date: 2012/05/15$
% $Revision: 1.1$  $Date: 2016/09/06$
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
disp('WSAF_DX_compare_demo');
% -------------------------------------------------------------------------

%%  Parameters setting

% Input parameters --------------------------------------------------------
Lx = 50000;                   % Length of the input signal
nRun = 10;                    % Number of runs
out_noise_level_dB = 60;      % SNR
out_noise_level = 10^(-out_noise_level_dB/20);   % Noise level

x = randn(Lx,1);              % x (Nx1) input signal array definition

% Colored signal generation -----------------------------------------------
a = 0.0;
b = sqrt(1-a^2);
%x = filter( b, [1 -a], randn( size(x))) ;  % H(z) = b/(1+a*z^-1)
disp('...... ');

% Adaptive filter definition ----------------------------------------------
M  =  10;             % Length of the linear filter
K  =   4;             % Affine Projection order; K=1 => NLMS
mu1 = 0.02;            % Learning rate of the linear filter
mQ1 = 0.02;            % Learning rate for ctrl points
delta = 0.1;          % APA regularization parameters
if Lx < 30000,        % Batch for evaluating MSE
        B = 100;
else
        B = 4000;
end

% Spline activation function definition and initialization ----------------
afinit  = 0;      % Init act. func. -1 0 ... (ONLY -1, bip.sig. or  0 =linear ) 
aftype  = 5;      % Kind act. f. -1 0 1 2 4 5; (4 = CR-spline, 5 = B-spline )
Slope   = 1;      % Slope
x_range = 2;      % Range limit
% Definition of several values of Delta X parameter -----------------------
DeltaX1 = 0.15;
DeltaX2 = 0.3;
DeltaX3 = 0.45;
DeltaX4 = 0.6;

% Creating the nonlinearity -----------------------------------------------
af1 = create_activation_function(afinit, aftype, DeltaX1, x_range, Slope, M);   % Delta X 1
af2 = create_activation_function(afinit, aftype, DeltaX2, x_range, Slope, M);   % Delta X 2
af3 = create_activation_function(afinit, aftype, DeltaX3, x_range, Slope, M);   % Delta X 3
af4 = create_activation_function(afinit, aftype, DeltaX4, x_range, Slope, M);   % Delta X 4


%% Initialization

% --- SAF definition ------------------------------------------------------
F1 = create_SPL_apa_adaptive_filter_1(M,K,mu1,mQ1,delta,af1);   % Delta X 1
F2 = create_SPL_apa_adaptive_filter_1(M,K,mu1,mQ1,delta,af2);   % Delta X 2
F3 = create_SPL_apa_adaptive_filter_1(M,K,mu1,mQ1,delta,af3);   % Delta X 3
F4 = create_SPL_apa_adaptive_filter_1(M,K,mu1,mQ1,delta,af4);   % Delta X 4


% Initialize --------------------------------------------------------------
N = Lx + M + 1;   % Total samples
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
disp('Algorithm start ... ');
t = clock;

for n = 0 : nRun-1
    fprintf('Test nr. %d/%d\n', n+1, nRun);
    
    % Nonlinear system signal generation ----------------------------------
    x = randn(Lx,1);           % x (Nx1) input signal array definition
    if a>0,
        x = filter( b, [1 -a], x);               % Colored input
    end
    z = filter([0.0154 0.0462 0.0465 0.0154], [1 -1.99 1.572 -0.4583], x);
    d = sin(z);
    %z = 0.1*(x.^3 - 3*x.^2 + x);                 % Nonlinearity
    %d = filter([1 0.75 0.5], [1 -0.75 0.5], z);  % Desired signal
    dn = out_noise_level*randn( size(x) );       % Noisy desired signal

    % SAF I.C. ------------------------------------------------------------
    F1.w(:) = 0;
    F1.w(1) = 1;  % Set filter i.c.
    F2.w(:) = 0;
    F2.w(1) = 1;  % Set filter i.c.
    F3.w(:) = 0;
    F3.w(1) = 1;  % Set filter i.c.
    F4.w(:) = 0;
    F4.w(1) = 1;  % Set filter i.c.
    
    % Set Activation Func I.C. --------------------------------------------
    F1.af.Q = af1.Q;
    F2.af.Q = af2.Q;
    F3.af.Q = af3.Q;
    F4.af.Q = af4.Q;
    
    % WSAF Evaluation -----------------------------------------------------  
    for k = 1 : Lx
        % Delta X 1 -------------------------------------------------------
        [F1, y1(k), e1(k)] = AF_APA_WSPL_F(F1, x(k), d(k) + dn(k) );
        
        % Delta X 2 -------------------------------------------------------
        [F2, y2(k), e2(k)] = AF_APA_WSPL_F(F2, x(k), d(k) + dn(k) );
        
        % Delta X 3 -------------------------------------------------------
        [F3, y3(k), e3(k)] = AF_APA_WSPL_F(F3, x(k), d(k) + dn(k) );
        
        % Delta X 4 -------------------------------------------------------
        [F4, y4(k), e4(k)] = AF_APA_WSPL_F(F4, x(k), d(k) + dn(k) );
    end
    
    % RUN time averaging
    em1 = em1 + e1.^2;    % Squared error
    em2 = em2 + e2.^2;    % Squared error
    em3 = em3 + e3.^2;    % Squared error
    em4 = em4 + e4.^2;    % Squared error
       
end

em1 = em1/nRun;    % MSE
em2 = em2/nRun;    % MSE
em3 = em3/nRun;    % MSE
em4 = em4/nRun;    % MSE
       
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
fprintf('HSAF with Delta X = %3.2f Steady-state MSE = %5.7f, equal to %5.3f dB\n',DeltaX1,mse1,10*log10(mse1));
fprintf('HSAF with Delta X = %3.2f Steady-state MSE = %5.7f, equal to %5.3f dB\n',DeltaX2,mse2,10*log10(mse2));
fprintf('HSAF with Delta X = %3.2f Steady-state MSE = %5.7f, equal to %5.3f dB\n',DeltaX3,mse3,10*log10(mse3));
fprintf('HSAF with Delta X = %3.2f Steady-state MSE = %5.7f, equal to %5.3f dB\n',DeltaX4,mse4,10*log10(mse4));
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
xa1 = zeros(1,KK);
dx = 2*xLIM/KK;
xx = -xLIM;   
for k = 1:KK
    yy1(k) = ActFunc(xx, F1.af);   % Delta X 1: Blue   
    yy2(k) = ActFunc(xx, F2.af);   % Delta X 2: Black
    yy3(k) = ActFunc(xx, F3.af);   % Delta X 3: Green
    yy4(k) = ActFunc(xx, F4.af);   % Delta X 4: Magenta
    xa1(k) = xx;
    xx = xx + dx;
end  
title('Profile of nonlinearities after learning','FontSize', 12, 'FontWeight', 'demi');
xlabel('Nonlinearity input {\itx}[{\itn}] ','FontSize', 12, 'FontWeight', 'demi');
ylabel('Nonlinearity output {\its}[{\itn}]','FontSize', 12, 'FontWeight', 'demi');
plot(xa1,yy1,'b','LineWidth',2);   % Delta X 1: Blue 
plot(xa1,yy2,'k','LineWidth',2);   % Delta X 2: Black
plot(xa1,yy3,'g','LineWidth',2);   % Delta X 3: Green
plot(xa1,yy4,'m','LineWidth',2);   % Delta X 4: Magenta
legend('Delta X 1','Delta X 2','Delta X 3','Delta X 4','Location','SouthEast');
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
plot(filter(bb,aa,edb1),'b','LineWidth',2);       % Delta X 1: Blue
plot(filter(bb,aa,edb2),'k','LineWidth',2);       % Delta X 2: Black
plot(filter(bb,aa,edb3),'g','LineWidth',2);       % Delta X 3: Green
plot(filter(bb,aa,edb4),'m','LineWidth',2);       % Delta X 4: Magenta
title('Comparisons of MSE by using different \Delta_x values','FontSize', 12, 'FontWeight', 'demi');
xlabel('Samples','FontSize', 12, 'FontWeight', 'demi');
ylabel('MSE [dB]   10log({\itJ}({\bfw},{\bfQ}))','FontSize', 12, 'FontWeight', 'demi');
legend('Delta X 1','Delta X 2','Delta X 3','Delta X 4');
set(gca, 'FontSize', 10, 'FontWeight', 'demi');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% © 2012
% M. Scarpiniti, D. Comminiello, R, Parisi and A. Uncini,  
% Department of Information Engineering, Electronics and Telecommunications
% (DIET) -- 'Sapienza' University of Rome
% -------------------------------------------------------------------------