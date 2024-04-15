function [coeffs, recreation] = Rectangle_Zernike_Fit(im, numZernike_start, numZernike_end)
% Rectangle_Zernike_Fit 对给定图像进行 Zernike 多项式拟合
%   im: 输入的图像矩阵
%   numZernike_start: 使用的 Zernike 多项式的起始序号
%   numZernike_end: 使用的 Zernike 多项式的结束序号
%   coeffs: 拟合得到的 Zernike 系数
%   recreation: 使用 Zernike 系数重建的图像

    % 确定单位圆内接矩形的位置
    [NumRows, NumCols] = size(im);
    k = NumCols / NumRows;
    if k >= 1
        x_max = sqrt(1 / (k^2 + 1));
        y_max = k * sqrt(1 / (k^2 + 1));
    else
        x_max = k * sqrt(1 / (k^2 + 1));
        y_max = sqrt(1 / (k^2 + 1));
    end

    % 生成笛卡尔坐标下的 Zernike 多项式
    [zern_expr, ~] = construct_zernike(numZernike_end);
    zern_xy = convertZernikePolarToCartesian(zern_expr);

    % 对矩形域下的 Zernike 多项式进行施密特正交化
    rect_zern = construct_rect_zernike(zern_xy, x_max, y_max);

    % 构造用于拟合的矩形域 Zernike 矩阵
    x = linspace(-x_max, x_max, NumCols);
    y = linspace(-y_max, y_max, NumRows);
    [x, y] = meshgrid(x, y);
    zernikeMatrices_xy = [];
    for i = numZernike_start:numZernike_end
        zernikeFunc = matlabFunction(rect_zern(i), 'Vars', {'x', 'y'});
        result = arrayfun(zernikeFunc, x, y);
        zernikeMatrices_xy(:,:,i) = reshape(result, size(x));
    end

    % 拟合 Zernike 系数
    z_mats_reshaped = [];
    for i = 1:size(zernikeMatrices_xy, 3)
        z_mats_reshaped(:, i) = reshape(zernikeMatrices_xy(:,:,i), NumCols * NumRows, 1);
    end
    im(isnan(im)) = 0;   
    z_mats_reshaped(isnan(z_mats_reshaped)) = 0;
    image_reshaped = reshape(im, NumCols * NumRows, 1);
    coeffs = (z_mats_reshaped.' * z_mats_reshaped) \ (z_mats_reshaped.' * image_reshaped);

    % 重建图像
    recreation = zeros(size(im));
    for i = 1:size(zernikeMatrices_xy, 3)
        recreation = recreation + coeffs(i) .* zernikeMatrices_xy(:,:,i);
    end
end







