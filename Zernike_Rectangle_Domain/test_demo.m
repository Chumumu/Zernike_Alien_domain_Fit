% 测试 Rectangle_Zernike_Fit 函数
clear
clc
clsoe all


% 步骤 1: 生成测试图像
% 生成笛卡尔坐标下的 Zernike 多项式
numZernike = 10; % 例如，使用 10 个 Zernike 多项式
[zern_expr, ~] = construct_zernike(numZernike);
zern_xy = convertZernikePolarToCartesian(zern_expr);
% 构造用于Zernike 矩阵
x = linspace(-1, 1, 100);
y = linspace(-1, 1, 100);
[x, y] = meshgrid(x, y);
zernikeMatrices_xy = [];
for i = 1:size(zern_xy,1)
    zernikeFunc = matlabFunction(zern_xy(i), 'Vars', {'x', 'y'});
    result = arrayfun(zernikeFunc, x, y);
    zernikeMatrices_xy(:,:,i) = reshape(result, size(x));
end
randomCoeffs = rand(1, numZernike) - 0.5;   % 生成随机 Zernike 系数
% 泽尼克合成波面
W = zeros(size(x)); % 初始化波前
for i = 1:numZernike
    W = W + randomCoeffs(i).*zernikeMatrices_xy(:,:,i);
end
% 裁剪出一个60*80的矩形区域用做拟合
Image = nan(size(x));
Image = W(20:80,10:90);



% 步骤 2: 调用 Rectangle_Zernike_Fit 函数
% 选择 Zernike 多项式的序号范围
numZernike_start = 1;
numZernike_end = 15;
[coeffs, recreation] = Rectangle_Zernike_Fit(Image, numZernike_start, numZernike_end);


% 步骤 3: 显示原始图像和重建图像
subplot(1, 3, 1);
imagesc(W); % 显示原始图像
colormap(jet); % 使用彩色图显示
colorbar;      % 显示颜色条
title('原始图像');
axis square;


subplot(1, 3, 2);
imagesc(Image); % 显示原始图像
colormap(jet); % 使用彩色图显示
colorbar;      % 显示颜色条
title('裁剪出来的待拟合面');


subplot(1, 3, 3);
imagesc(recreation); % 显示重建图像
colormap(jet); % 使用彩色图显示
colorbar;      % 显示颜色条
title('重建后的拟合面');



fitError = recreation - Image;
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