function Z_cartesian = convertZernikePolarToCartesian(Z_polar)
% convertZernikePolarToCartesian 将 Zernike 多项式从极坐标系转换到笛卡尔坐标系
%   Z_polar: 极坐标系下的 Zernike 多项式数组
%   Z_cartesian: 笛卡尔坐标系下的 Zernike 多项式数组

    syms r t x y
    Z_cartesian = sym(zeros(size(Z_polar))); % 初始化笛卡尔坐标系下的 Zernike 多项式向量

    for i = 1:length(Z_polar)
        Z_polar_expanded = expand(Z_polar(i)); % 展开每个多项式
        % 替换为笛卡尔坐标
        Z_cartesian(i) = simplify(subs(Z_polar_expanded, [r*cos(t), r*sin(t), r^2], [x, y, x^2 + y^2])); 
    end
end

