function [annulus_zern_expr] = construct_annulus_zernike(zern_expr, m)
%   zern_expr : 圆域zernike多项式
%   m : 圆环的内径
%   annulus_zern_expr ： zern_expr经过Gram-Schimdt正交化后得到的环域Zernike多项式'
%     syms m % 如果要看看多项式的正确性取消注释
    num = length(zern_expr);
    A = sym(zeros(num,num));   % A为系数矩阵
    A(1,1) = 1 ;
    
    % 1、正交化
    X = zern_expr(1);
    for i = 2: num
        T = sym(zeros(i-1,i-1));
        H = sym(zeros(i-1,1));
        for p = 1:i-1
            for q = 1:i-1
                T(p,q) = annulus_integral(X(p),zern_expr(q),m);
            end
            H(p) = -annulus_integral(X(p),zern_expr(i),m);
        end
        % 求解线性方程组 Tx = H
        x = T \ H;
        % 需要把系数矩阵存入A中
        x = sym([x; 1]);
        A(i,1:i) = x;
        X = A(1:i,1:i)*(zern_expr(1:i));
        X = simplify(X);
    end    
    
    % 2、单位化
    annulus_zern_expr = sym(zeros(size(X)));
    for i = 1:length(annulus_zern_expr)
        annulus_zern_expr(i) =  X(i)./simplify((annulus_integral(X(i),X(i),m)).^(1/2));
        annulus_zern_expr(i) = simplify(annulus_zern_expr(i));
    end
end

function [Q] = annulus_integral(f1,f2,m)
%   输入为函数a和函数b，输出为函数乘积f在环形域内的积分Q
syms r t 

    
    % 定义积分函数
    f = f1 .* f2.*r;
    s = pi*(1 - m^2 );
    % 执行双重积分
    F = int(int(f, r, m, 1), t, 0, 2*pi);
    F = F ./s;
    Q = simplify(F);
    % disp(Q);

end

