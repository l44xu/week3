clear;clc
% %% Direct inversion
G = 6.674e-11;
m1 = 10;
m2 = 50;
r = 100; 
F = G*m1*m2/r^2;

r_hat = (G*m1*m2/F)^0.5;

r - r_hat

%% Look up table (LUT) approach
% define the LUT
[m1,m2,r] = ndgrid(1:100,1:100,1:10);
m1 = m1(:);
m2 = m2(:);
r = r(:);
F = (G*m1.*m2)./(r.^2);
x = [m1,m2,F];
% define new observations
m1_test = [rand(50,1)*100; rand(50,1)*100+100];
m2_test = [rand(50,1)*100; rand(50,1)*100+100];
r_true = [rand(50,1)*10; rand(50,1)*100+10];
F_test = (G*m1_test.*m2_test)./(r_true.^2);

test.m1=m1_test;
test.m2=m2_test;
test.r=r_true;
test.F = F_test;
save('test','test');


E_error = zeros(100,1);
M_error = zeros(100,1);
E_r_hat = zeros(100,1);
M_r_hat = zeros(100,1);
var_x = var(x);
for i = 1 : 100
    x_test_i = [m1_test(i),m2_test(i),F_test(i)];
    % metric 1 : Euclidean distance
    E_dists = sum((x-repmat(x_test_i,100^2*10,1)).^2, 2); % Euclidean distance
    [~,min_idx_E] = min(E_dists); % find the most similar combination
    E_r_hat(i) = r(min_idx_E);
    E_error(i) = (r_true(i) - r(min_idx_E))^2;
    
    % metric 2 : Mahalanobis distance
    M_dists = sum((x-repmat(x_test_i,100^2*10,1)).^2./(repmat(var_x,100^2*10,1)), 2); % M distance
    [~,min_idx_M] = min(M_dists); % find the most similar combination
    M_error(i) = (r_true(i) - r(min_idx_M))^2;
    M_r_hat(i) = r(min_idx_M);
end

figure; plot([E_error(:),M_error(:)]); legend('E distance', 'M distance');
figure; plot([r_true, E_r_hat]); legend('r true', 'r E estimate'); title('LUT E results'); 
figure; plot([r_true, M_r_hat]); legend('r true', 'r M estimate'); title('LUT M results'); 

% (1) interpolation good!
% (2) extrapolation bad!
% (3) Mahalanobis better than Euclidean
% (4) other metrics
% (5) choose more candidates and use their mean or median

%% Artificial Neural Network (ANN) trained on simulated data from physical models
% generate simulated data
[m1,m2,r] = ndgrid(1:100,1:100,1:10);
m1 = m1(:)';
m2 = m2(:)';
r = r(:)';
F = (G*m1.*m2)./(r.^2);

inputs = [F;m1;m2];
targets = r;

% Create a Fitting Network
hiddenLayerSize = 2;
net = fitnet(hiddenLayerSize,'trainlm');
% net.layers{1}.transferFcn = 'logsig';
% net.layers{2}.transferFcn = 'logsig';

% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
 
% Train the Network
[net,tr] = train(net,inputs,targets);
 
% Test the Network
outputs = net(inputs);
errors = gsubtract(outputs,targets);
performance = perform(net,targets,outputs);
 
% View the Network
view(net)

%% interpolation & extrapolation
load test;
rhat = net([test.F';test.m1';test.m2']);
figure; plot([test.r, rhat']); legend('r true', 'r estimate'); title('ANN one layer (2 neurons) results'); 

%% Artificial Neural Network (ANN) trained on simulated data from physical models
% generate simulated data
[m1,m2,r] = ndgrid(1:100,1:100,1:10);
m1 = m1(:)';
m2 = m2(:)';
r = r(:)';
F = (G*m1.*m2)./(r.^2);

inputs = [F;m1;m2];
targets = r;

% Create a Fitting Network
hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize,'trainlm');
% net.layers{1}.transferFcn = 'logsig';
% net.layers{2}.transferFcn = 'logsig';

% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
 
% Train the Network
[net,tr] = train(net,inputs,targets);
 
% Test the Network
outputs = net(inputs);
errors = gsubtract(outputs,targets);
performance = perform(net,targets,outputs);
 
% View the Network
view(net)

%% interpolation & extrapolation
load test;
rhat = net([test.F';test.m1';test.m2']);
figure; plot([test.r, rhat']); legend('r true', 'r estimate'); title('ANN one layer (10 neurons) results'); 


%% MLE
% get the observed data
m1 = 100 + randn(100,1);
m2 = 200 + randn(100,1)*2;
r = 50;
F = (G*m1.*m2)./(r^2);
% estimate using MLE
r_hat = ((sum(G*m1.^2.*m2.^2))/(sum(F.*m1.*m2)))^0.5;
error = r - r_hat;

%% MAP
% case 1: 'r' is model parameter, which is a random variable satisfying uniform
% distribution from 'min_r' to 'max_r'
fun = @(r) sum((F-G*m1.*m2/r^2).^2);
min_r = 20;
max_r = 70;
r0 = 65;
r_hat = simulannealbnd(fun,r0,min_r,max_r);

% case 2: 'r' is model parameter, which is a random variable satisfying Gaussian
% distribution with 'mu=50', 'sigma=2';
mu_r = 50;
sigma_r = 2;
r0 = 80;
beta = 0.1;
fun = @(r) sum((F-G*m1.*m2/r^2).^2)+beta*(r-mu_r)^2/(sigma_r^2);
r_hat = simulannealbnd(fun,r0);

% case 3: G is model parameter, 'r' is laten variable, which is a random variable 
% satisfying a mixture of two Gaussian distributions with known model
% parameters
% step (1): simulate data according forward model
mu1 = 20;
mu2 = 60;
sigma1 = 0.1;
sigma2 = 0.2;
p1 = 0.5;
p2 = 0.5;

m1 = 100 + randn(100,1);
m2 = 200 + randn(100,1)*2;
r = [randn(50,1)*sigma1+mu1;randn(50,1)*sigma2+mu2];
r = r(randperm(length(r)));
F = (G*m1.*m2)./(r.^2);

% step (2): perform data inversion with MAP 
% inite estimate of 'r'
r_hat = rand(100,1)*100;
l(r_hat>50) = 2;
l(r_hat<=50) = 1;
theta = zeros(100,5);
for i = 1:100
    i
    % M-step: estimate 'G' using r based on MLE
    G_hat = sum((F.*m1.*m2)./r_hat.^2)/sum((m1.^2.*m2.^2)./r_hat.^4); theta(i,1) = G_hat;
    mu1_hat = mean(r_hat(l==1));  theta(i,2) = mu1_hat;
    sigma1_hat = std(r_hat(l==1)); theta(i,3) = sigma1_hat;
    mu2_hat = mean(r_hat(l==2)); theta(i,4) = mu2_hat;
    sigma2_hat = std(r_hat(l==2)); theta(i,5) = sigma2_hat;
    theta(i,:)
    G_hat
    % E-step: estimate 'r' using 'G'
    % (1) estimate l given r
    dist1 = ((r_hat-mu1_hat).^2/sigma1_hat^2);
    dist2 = ((r_hat-mu2_hat).^2/sigma2_hat^2);
    dists = [dist1,dist2];
    [~,l] = min(dists,[],2);
    % (2) estimate r given l
    beta = 10;
    fun = @(r_hat) sum((F-G_hat*m1.*m2./r_hat.^2).^2)+beta*sum(((r_hat-mu1_hat).^2/sigma1_hat^2).*(l==1) + ((r_hat-mu2_hat).^2/sigma2_hat^2).*(l==2));
    r_hat = simulannealbnd(fun,r_hat);
end

figure; plot(theta); legend('G','mu1','sigma1','mu2','sigma2');

dist1 = ((r-mu1).^2/sigma1^2);
dist2 = ((r-mu2).^2/sigma2^2);
dists = [dist1,dist2];
[~,l] = min(dists,[],2);

dist1 = ((r_hat-mu1).^2/sigma1^2);
dist2 = ((r_hat-mu2).^2/sigma2^2);
dists = [dist1,dist2];
[~,l_hat] = min(dists,[],2);

figure;imagesc(reshape(l,10,10));
figure;imagesc(reshape(l_hat,10,10));



