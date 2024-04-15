function [zern_expr,indices] = construct_zernike(num)
% 构造圆域zernike多项式的前num项
syms r t
zern_expr = [];
indices = [];

for n = 0:99
    for m = -n:2:n
        zern_expr = [zern_expr;simplify(zernike(r,t,n,m))];   % Zernike表达式
        indices = [indices; n m];                             % Zernike表达式对应的m，n索引
        if size(indices,1) == num
          return
        end
    end
end
end


function zern = zernike(r,t,n,m)
% radial degree       : n
% azimuthal frequency : m 
% 对应的Zernike多项式
    if mod(n-m,2) == 1
        error('n-m must be even');
    end
    if n < 0
        error('n must both be positive')
    end
    if floor(n) ~= n || floor(m) ~= m
        error('n and m must both be integers')
    end
    

    % (这里加上对Zernike多项式进行归一化）
    if m < 0
        zern = -(2*n+2)^(1/2).*zernike_radial(r,n,-m).*sin(m*t);
    elseif m > 0
        zern = (2*n+2)^(1/2).*zernike_radial(r,n,m).*cos(m*t);
    else
        zern = (n+1).^(1/2).*zernike_radial(r,n,m);
    end
end


function radial = zernike_radial(r,n,m)
% Zernike多项式的径向部分
    if mod(n-m,2) == 1
        error('n-m must be even');
    end
    if n < 0 || m < 0
        error('n and m must both be positive in radial function')
    end
    if floor(n) ~= n || floor(m) ~= m
        error('n and m must both be integers')
    end
    if n == m
        radial = r.^n;
    elseif n - m == 2
        radial = n*zernike_radial(r,n,n)-(n-1)*zernike_radial(r,n-2,n-2);
    else
        H3 = (-4*((m+4)-2)*((m+4)-3)) / ((n+(m+4)-2)*(n-(m+4)+4));
        H2 = (H3*(n+(m+4))*(n-(m+4)+2)) / (4*((m+4)-1))  +  ((m+4)-2);
        H1 = ((m+4)*((m+4)-1) / 2)  -  (m+4)*H2  +  (H3*(n+(m+4)+2)*(n-(m+4))) / (8);
        radial = H1*zernike_radial(r,n,m+4) + (H2+H3 ./ r.^2).*zernike_radial(r,n,m+2);
    end

end