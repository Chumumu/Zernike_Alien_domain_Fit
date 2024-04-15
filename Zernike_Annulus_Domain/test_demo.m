% 测试 Rectangle_Zernike_Fit 函数
clear
clc
close all


% 步骤 1: 生成测试图像
% 生成笛卡尔坐标下的 Zernike 多项式
m = 0.5; % 内径设为0.5
numZernike = 10; % 例如，使用 10 个 Zernike 多项式
[zern_expr, ~] = construct_zernike(numZernike);

% 构造用于Zernike 矩阵
x = linspace(-1, 1, 100);
y = linspace(-1, 1, 100);
[x, y] = meshgrid(x, y);
[t, r] = cart2pol(x, y);
zernikeMatrices = [];
for i = 1:size(zern_expr,1)
    zernikeFunc = matlabFunction(zern_expr(i), 'Vars', {'r', 't'});
    result = arrayfun(zernikeFunc, r, t);
    zernikeMatrices(:,:,i) = reshape(result, size(x));
    zernikeMatrices(:,:,i) = crop_annulus(zernikeMatrices(:,:,i), m,1);
end
randomCoeffs = rand(1, numZernike) - 0.5;   % 生成随机 Zernike 系数
% 泽尼克合成波面
W = zeros(size(x)); % 初始化波前
for i = 1:numZernike
    W = W + randomCoeffs(i).*zernikeMatrices(:,:,i);
end




% 步骤 2: 调用 Rectangle_Zernike_Fit 函数
% 选择 Zernike 多项式的序号范围
numZernike_start = 1;
numZernike_end = 10;
[coeffs, recreation] = Annulus_Zernike_Fit(W,m, numZernike_start, numZernike_end);


% 步骤 3: 显示原始图像和重建图像
subplot(1, 2, 1);
imagesc(W); % 显示原始图像
colormap(jet); % 使用彩色图显示
colorbar;      % 显示颜色条
title('原始图像');
axis square;

subplot(1, 2, 2);
imagesc(recreation); % 显示重建图像
colormap(jet); % 使用彩色图显示
colorbar;      % 显示颜色条
title('重建后的拟合面');



fitError = recreation - W;
RMS_daoru_gailiang = std(fitError,'omitnan');
PV_daoru_gailiang = max(max(fitError))-min(min(fitError));
figure;
surf(fitError,'EdgeColor','none')
set(gcf,'color','white');
daspect([1 1 1]);
grid off;
shading interp; %颜色渐变
h=colorbar;
% set(get(h,'Title'),'string','m'); % 此处根据需要设置colorbar的单位【可修改】
title('最小二乘法所求面形与导入原始面形的偏差');
view(2);
xlim([0 size(fitError,2)]);
ylim([0 size(fitError,1)]);