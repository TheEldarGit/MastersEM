%%%%%%%%%%%%%% Selected Topics in Image Processing %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PROJECT : Implementation of Chrochrome On ViaReggio Data Set %%%%%%
% Eldar Mamedov    I.D. 313251043
% Yaniv Lerner     I.D. 312218837
%% Clean project 
close all;clear all; clc; 
%% Redaing the files 
%  D1F12H1– D2F22H2 
% will read the 2 files to detect anumalys

day1_header_path = "D1_F12_H1_Cropped\Full_spectral_samples\Destriped_data\D1_F12_H1_Cropped_des.hdr";
day1_data_path = "D1_F12_H1_Cropped\Full_spectral_samples\Destriped_data\D1_F12_H1_Cropped_des";
day2_header_path = "D2_F22_H2_Cropped_Aligned\Full_spectral_samples\Destriped_data\D2_F22_H2_Cropped_des_Aligned.hdr"; 
day2_data_path = "D2_F22_H2_Cropped_Aligned\Full_spectral_samples\Destriped_data\D2_F22_H2_Cropped_des_Aligned";

day1_2_header_path = "D1_F12_H2_Cropped_Aligned\Full_spectral_samples\Destriped_data\D1_F12_H2_Cropped_des_Aligned.hdr"; 
day1_2_data_path = "D1_F12_H2_Cropped_Aligned\Full_spectral_samples\Destriped_data\D1_F12_H2_Cropped_des_Aligned";

path_graphes = "Images destripte\";

% day1_header_path = "D1_F12_H1_Cropped\Full_spectral_samples\Noise_whitened_data\D1_F12_H1_Cropped_NW.hdr";
% day1_data_path = "D1_F12_H1_Cropped\Full_spectral_samples\Noise_whitened_data\D1_F12_H1_Cropped_NW";
% day2_header_path = "D2_F22_H2_Cropped_Aligned\Full_spectral_samples\Noise_whitened_data\D2_F22_H2_Cropped_NW_Aligned.hdr"; 
% day2_data_path = "D2_F22_H2_Cropped_Aligned\Full_spectral_samples\Noise_whitened_data\D2_F22_H2_Cropped_NW_Aligned";
% 
% day1_2_header_path = "D1_F12_H2_Cropped_Aligned\Full_spectral_samples\Noise_whitened_data\D1_F12_H2_Cropped_NW_Aligned.hdr"; 
% day1_2_data_path = "D1_F12_H2_Cropped_Aligned\Full_spectral_samples\Noise_whitened_data\D1_F12_H2_Cropped_NW_Aligned";
%
% path_graphes = "Images white noise\";

% X is refarnce and Y is test 
hcube_X = hypercube(day1_data_path,day1_header_path);
hcube_Y = hypercube(day2_data_path,day2_header_path);
% X is refarnce and W is test
hcube_W = hypercube(day1_2_data_path,day1_2_header_path);

%% visualizing the data
X_cube = hcube_X.DataCube;
X_img = im2double(X_cube);
[rows_1,cols_1,labels_day1] = size(X_cube);

Y_cube = hcube_Y.DataCube;
Y_img = im2double(Y_cube);
[rows_2,cols_2,labels_day2] = size(Y_img);

W_cube = hcube_W.DataCube;
W_img = im2double(W_cube);
[rows_2,cols_2,labels_day2] = size(W_img);

%band = [22,32,64];

rgbImg_X =  colorize(hcube_X,"Method","rgb","ContrastStretching",true);
rgbImg_Y =  colorize(hcube_Y,"Method","rgb","ContrastStretching",true);
rgbImg_W =  colorize(hcube_W,"Method","rgb","ContrastStretching",true);
%%
figure; imagesc(rgbImg_X);title('RGB Image X');
figure; imagesc(rgbImg_Y);title('RGB Image Y');
figure; imagesc(rgbImg_W);title('RGB Image W');

%% mean and Covarince
X_reshaped_org = reshape(X_cube,rows_1*cols_1,labels_day1);
Y_reshaped_org = reshape(Y_cube,rows_1*cols_1,labels_day1);
W_reshaped_org = reshape(W_cube,rows_1*cols_1,labels_day1);

X_reshaped_withoutMean = X_reshaped_org - mean(X_reshaped_org);
Y_reshaped_withoutMean = Y_reshaped_org - mean(Y_reshaped_org);
W_reshaped_withoutMean = W_reshaped_org - mean(W_reshaped_org);

% Covarince 
cov_Mat_X = X_reshaped_withoutMean'*X_reshaped_withoutMean;
cov_Mat_X = cov_Mat_X./(rows_1*cols_1);

cov_Mat_XY = X_reshaped_withoutMean'*Y_reshaped_withoutMean;
cov_Mat_XY = cov_Mat_XY./(rows_1*cols_1);

cov_Mat_XW = X_reshaped_withoutMean'*W_reshaped_withoutMean;
cov_Mat_XW = cov_Mat_XW./(rows_1*cols_1);


cov_Mat_sub = (Y_reshaped_withoutMean-X_reshaped_withoutMean)'*(Y_reshaped_withoutMean-X_reshaped_withoutMean)./(rows_1*cols_1);
cov_Mat_sub_W = (W_reshaped_withoutMean-X_reshaped_withoutMean)'*(W_reshaped_org-X_reshaped_withoutMean)./(rows_1*cols_1);

%% Create L
L = cov_Mat_XY*inv(cov_Mat_X); 
L_W = cov_Mat_XW*inv(cov_Mat_X); 
%% y -Lx
e_reshaped = (Y_reshaped_withoutMean' - L * X_reshaped_withoutMean')';
E_cov = cov(e_reshaped);
e = reshape(e_reshaped, rows_1,cols_1,labels_day1);

e_reshaped_W = (W_reshaped_withoutMean' - L_W * X_reshaped_withoutMean')';
E_cov_W = cov(e_reshaped_W);
e_W = reshape(e_reshaped_W, rows_1,cols_1,labels_day1);


e_sub = (Y_reshaped_withoutMean - X_reshaped_withoutMean);
e_sub_w = (W_reshaped_withoutMean - X_reshaped_withoutMean);

%% anomaly rate 
A_rate = zeros(1,rows_1*cols_1);
A_rate_W = zeros(1,rows_1*cols_1);
A_rate_sub =  zeros(1,rows_1*cols_1);
A_rate_sub_w =  zeros(1,rows_1*cols_1);
E_cov_inv = E_cov^(-1);
E_cov_inv_W = E_cov_W^(-1);
E_cov_sub_inv = cov_Mat_sub^(-1);
E_cov_sub_inv_W= cov_Mat_sub_W^(-1);

for jj = 1:rows_1*cols_1
    A_rate(jj) = e_reshaped(jj,:)*E_cov_inv*(e_reshaped(jj,:))';
    A_rate_W(jj) = e_reshaped_W(jj,:)*E_cov_inv_W*(e_reshaped_W(jj,:))';
    A_rate_sub(jj) = e_sub(jj,:)*E_cov_sub_inv*e_sub(jj,:)';
    A_rate_sub_W(jj) = e_sub_w(jj,:)*E_cov_sub_inv_W*e_sub_w(jj,:)';
end

mean(A_rate,"all")
mean(A_rate_W,"all")
mean(A_rate_sub,"all")
mean(A_rate_sub_W,"all")
%% Loading true Mask
% load('CHANGE_REFERENCE_MAPS\TEST_D2_F22_H2_REF_D1_F12_H1\MASK_TEST_D2_F22_H2_REF_D1_F12_H1.mat');
% load('CHANGE_REFERENCE_MAPS\TEST_D1_F12_H1_REF_D2_F22_H2\MASK_TEST_D1_F12_H1_REF_D2_F22_H2.mat');
% figure;
% imshow(MASK_TEST_D2H2_REF_D1H1');
% figure;
% imshow(MASK_TEST_D1H1_REF_D2H2');

%% 
A_rate_reshape = reshape(A_rate, rows_1,cols_1);
A_rate_reshape_W = reshape(A_rate_W, rows_1,cols_1);
A_rate_sub_reshape = reshape(A_rate_sub, rows_1,cols_1);
A_rate_sub_reshape_W = reshape(A_rate_sub_W, rows_1,cols_1);
mult_parm = 3;
thrashold = mean(A_rate_reshape,"all")+mult_parm*std(A_rate_reshape(:));
thrashold_W = mean(A_rate_reshape_W,"all")+mult_parm*std(A_rate_reshape_W(:));
thrashold_sub = mean(A_rate_sub_reshape,"all")+mult_parm*std(A_rate_sub_reshape(:));
thrashold_sub_w = mean(A_rate_sub_reshape_W,"all")+mult_parm*std(A_rate_sub_reshape_W(:));
% Using the algorithem 
RGB_1 = insertObjectMask(rgbImg_X,abs(A_rate_reshape)>thrashold);

% Groud true 
RGB_2 = insertObjectMask(rgbImg_X,abs(A_rate_reshape_W)> thrashold_W);
RGB_3 = insertObjectMask(rgbImg_X,abs(A_rate_sub_reshape)>thrashold_sub);
RGB_4 = insertObjectMask(rgbImg_X,abs(A_rate_sub_reshape_W)>thrashold_sub_w);

figure;
subplot(1,2,1)
imshow(RGB_1); title("X with chrochrome map day2")
subplot(1,2,2)
imshow(RGB_3);title("X with map of Y -X ")
figure;
subplot(1,2,1)
imshow(RGB_2); title("X with chrochrome map day1")
subplot(1,2,2)
imshow(RGB_4); title("X with sub map day1")
%% histogram A rate 
figure;
histogram(A_rate_reshape); title("Anomaly X-> Yrate histogram")
xline(thrashold,'-',{'mean +3*std'})
figure;
histogram(A_rate_sub); title("Anomaly X->Y sub rate histogram")
xline(thrashold_sub,'-',{'mean +3*std'})

figure;
histogram(A_rate_reshape_W); title("Anomaly X->W rate histogram")
xline(thrashold_W,'-',{'mean +3*std'})

%% 
A_rate_reshape = reshape(A_rate, rows_1,cols_1);
title_save = "Anomaly rate";
figure;
imshow((A_rate_reshape-min(A_rate_reshape(:)))/(max(A_rate_reshape(:))-min(A_rate_reshape(:))));

title("Anomaly rate");
colormap('jet')
colorbar;
%% inplanting target and creating roc curve 
Y_org= reshape(Y_reshaped_org,rows_1,cols_1,labels_day1);
t = squeeze(Y_org(25,25,:));


W_org= reshape(W_reshaped_org,rows_1,cols_1,labels_day1);
t_W = squeeze(W_org(25,25,:));

p = 0.2;

Y_WT = Y_reshaped_withoutMean + p*t';
W_WT = W_reshaped_withoutMean + p*t_W';

%%

e_reshaped_tag = ((Y_WT)' - L * (X_reshaped_withoutMean)')';
e_reshaped_tag_sub = (Y_WT -X_reshaped_withoutMean);
e_reshaped_tag_sub_W = (W_WT -X_reshaped_withoutMean);
e_reshaped_W_tag = ((W_WT)' - L_W * (X_reshaped_withoutMean)')';
% anomaly rate with target 
A_rate_tag = zeros(1,rows_1*cols_1);
A_rate_tag_sub = zeros(1,rows_1*cols_1);
A_rate_tag_W = zeros(1,rows_1*cols_1);
A_rate_tag_sub_W = zeros(1,rows_1*cols_1);
for jj = 1:rows_1*cols_1
    A_rate_tag(jj) = e_reshaped_tag(jj,:)*E_cov_inv*(e_reshaped_tag(jj,:))';
    A_rate_tag_W(jj) = e_reshaped_W_tag(jj,:)*E_cov_inv_W*(e_reshaped_W_tag(jj,:))';
    A_rate_tag_sub(jj) = e_reshaped_tag_sub(jj,:)*E_cov_sub_inv*(e_reshaped_tag_sub(jj,:))';
    A_rate_tag_sub_W(jj) = e_reshaped_tag_sub_W(jj,:)*E_cov_sub_inv_W*(e_reshaped_tag_sub_W(jj,:))';
end
mean(A_rate_tag,"all")
mean(A_rate_tag_W,"all")
mean(A_rate_tag_sub,"all")


%%  Histograms

[NT_val,NT_bins] = histcounts(A_rate,1000);
[WT_val,WT_bins] = histcounts(A_rate_tag,1000);

figure;
plot(NT_bins(1:1000),NT_val);

hold on;
plot(WT_bins(1:1000),WT_val);
title('Original Histogram');
legend('No target','With target');
hold off;

%% Roc Curve

[X,Y,T,AUC] = ROC_EM(A_rate,A_rate_tag);
[W,Z,T_sub,AUC_sub] = ROC_EM(A_rate_sub,A_rate_tag_sub);
[R,V,T_PCA,AUC_day1] = ROC_EM(A_rate_W,A_rate_tag_W);
[X,Q,T_PCA,AUC_day1_sub] = ROC_EM(A_rate_sub_W,A_rate_tag_sub_W);
figure; 
plot(X,Y); title("ROC: Chrochrome-Subtraction")
grid on
hold on
plot(W,Z)
plot(R,V)
plot(X,Q)
% legend(name)
xlabel('Pfa');
ylabel('Pd');
legend('Chrochrome','SUB','Chrochrome day1','SUB day1')

hold off
