function rect_zern_expr = construct_rect_zernike(Zernike_xy, x_limit, y_limit)
% construct_rect_zernike 在矩形域内对笛卡尔坐标系下的 Zernike 多项式进行施密特正交化
%   Zernike_xy: 圆域 Zernike 多项式表达式数组
%   x_limit, y_limit: 矩形域二重积分 x 和 y 的上下限
%   rect_zern_expr: 施密特正交化后得到的矩形域 Zernike 多项式数组

    syms x y real
    num = length(Zernike_xy);
    A = sym(zeros(num, num)); % 系数矩阵
    A(1, 1) = 1;

    % 正交化
    X = Zernike_xy(1);
    for i = 2:num
        T = sym(zeros(i-1, i-1));
        H = sym(zeros(i-1, 1));
        for p = 1:i-1
            for q = 1:i-1
                T(p, q) = rect_integral(X(p), Zernike_xy(q), x_limit, y_limit);
            end
            H(p) = -rect_integral(X(p), Zernike_xy(i), x_limit, y_limit);
        end
        % 求解线性方程组 Tx = H
        x = T \ H;
        % 存入系数矩阵
        x = sym([x; 1]);
        A(i, 1:i) = x;
        X = A(1:i, 1:i) * Zernike_xy(1:i);
        X = simplify(X);
    end    
    % 单位化
    rect_zern_expr = sym(zeros(size(X)));
    for i = 1:length(rect_zern_expr)
        rect_zern_expr(i) = X(i) ./ simplify((rect_integral(X(i), X(i), x_limit, y_limit))^(1/2));
        rect_zern_expr(i) = simplify(rect_zern_expr(i));
    end
end

function Q = rect_integral(f1, f2, x_limit, y_limit)
% rect_integral 计算矩形域内两个函数的内积
%   f1, f2: 两个函数表达式
%   x_limit, y_limit: 矩形域二重积分 x 和 y 的上下限
%   Q: 函数 f1 和 f2 在矩形域内的内积

    syms x y real
    % 定义乘积函数
    f = f1 * f2;
    % 执行双重积分
    s = 4 * x_limit * y_limit;
    F = int(int(f, x, -x_limit, x_limit), y, -y_limit, y_limit);
    F = F / s;
    Q = simplify(F);
end
