function [coeffs, recreation] = Annulus_Zernike_Fit(im, m, numZernike_start, numZernike_end)
    % Annulus_Zernike_Fit 对给定图像进行环域 Zernike 多项式拟合
    % 参数:
    %   im: 输入的图像矩阵
    %   m: 遮光比（也即内径）
    %   numZernike_start: 使用的 Zernike 多项式的起始序号
    %   numZernike_end: 使用的 Zernike 多项式的结束序号
    % 返回值:
    %   coeffs: 拟合得到的 Zernike 系数
    %   recreation: 使用 Zernike 系数重建的图像

    % 检查输入图像是否为方阵
    [NumRows, NumCols] = size(im);
    if NumRows ~= NumCols
        error('输入图像行数和列数必须相等');
    end

    % 裁剪输入的图像为环域
    x = linspace(-1, 1, NumCols);
    y = linspace(-1, 1, NumRows);
    [x, y] = meshgrid(x, y);
    [t, r] = cart2pol(x, y);
    im = crop_annulus(im, m, 1);  % 使用 crop_annulus 裁剪
    
    % 生成 Zernike 多项式
    [zern_expr, ~] = construct_zernike(numZernike_end);

    % 经过 Gram-Schmidt 正交化后得到的内径比为 m 的环域 Zernike 多项式
    annulus_zern = construct_annulus_zernike(zern_expr, m);

    % 构造用于拟合的环域 Zernike 矩阵
    zernikeMatrices = [];
    for i = numZernike_start:numZernike_end
        zernikeFunc = matlabFunction(annulus_zern(i), 'Vars', {'t', 'r'});
        result = arrayfun(zernikeFunc, t, r);
        zernikeMatrices(:, :, i) = reshape(result, size(x));
        zernikeMatrices(:, :, i) = crop_annulus(zernikeMatrices(:, :, i), m, 1);  % 对每项 Zernike 矩阵进行裁剪
    end

    % 拟合 Zernike 系数
    z_mats_reshaped = [];
    for i = 1:size(zernikeMatrices, 3)
        z_mats_reshaped(:, i) = reshape(zernikeMatrices(:, :, i), NumCols * NumRows, 1);
    end
    im(isnan(im)) = 0;
    z_mats_reshaped(isnan(z_mats_reshaped)) = 0;
    image_reshaped = reshape(im, NumCols * NumRows, 1);
    coeffs = (z_mats_reshaped.' * z_mats_reshaped) \ (z_mats_reshaped.' * image_reshaped);

    % 重建图像
    recreation = zeros(size(im));
    for i = 1:size(zernikeMatrices, 3)
        recreation = recreation + coeffs(i) * zernikeMatrices(:, :, i);
    end
end







