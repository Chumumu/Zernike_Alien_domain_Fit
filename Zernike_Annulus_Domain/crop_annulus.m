function cropped_im = crop_annulus(im, inner_frac, outer_frac)
    % crop_annulus 裁剪图像为环域形状
    % 参数:
    %   im: 被裁剪图像
    %   inner_frac: 内径
    %   outer_frac: 外径（默认为1）
    % 返回值:
    %   cropped_im: 裁剪后的图像
    
    if nargin < 3
        outer_frac = 1;
    end

    if inner_frac < 0 || inner_frac > 1 || outer_frac < 0 || outer_frac > 1
        error('inner_frac and outer_frac must have values between 0 and 1');
    end
    
    cropped_im = im;
    center_x = (size(im, 2) + 1) / 2;
    center_y = (size(im, 1) + 1) / 2;
    inner_radius = (size(im, 2) - center_x) * inner_frac;
    outer_radius = (size(im, 2) - center_x) * outer_frac;

    for row = 1:size(im, 1)
        for col = 1:size(im, 2)
            distance_to_center = sqrt((row - center_y)^2 + (col - center_x)^2);
            if distance_to_center < inner_radius || distance_to_center > outer_radius
                cropped_im(row, col) = nan;
            end
        end
    end
end