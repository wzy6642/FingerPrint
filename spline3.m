function [] = spline3( X, Y, dY, x0, m )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
N = size(X, 2);
s0 = dY(1);
sN = dY(2);
interval = 0.025;
disp('x0为插值点')
x0
h = zeros(1, N-1);
for i=1:N-1
    h(1,i)=X(i+1)-X(i);
end
d(1,1) = 6*((Y(1,2)-Y(1,1))/h(1,1)-s0)/h(1,1);
d(N,1) = 6*(sN-(Y(1,N)-Y(1,N-1))/h(1,N-1))/h(1,N-1);
for i=2:N-1
    d(i,1)=6*((Y(1,i+1)-Y(1,i))/h(1,i)-(Y(1,i)-Y(1,i-1))/h(1,i-1))/(h(1,i)+h(1,i-1));
end
mu = zeros(1,N-1);
md = zeros(1,N-1);
md(1,N-1) = 1;
mu(1,1) = 1;
for i = 1:N-2
    u = h(1,i+1)/(h(1,i)+h(1,i+1));
    mu(1,i+1) = u;
    md(1,i) = 1-u;
end
p(1,1) = 2;
q(1,1) = mu(1,1)/2;
for i = 2:N-1
    p(1,i) = 2 - md(1,i-1)*q(1,i-1);
    q(1,i) = mu(1,i) / p(1,i);
end
p(1,N) = 2-md(1,N-1)*q(1,N-1);
y = zeros(1,N);
y(1,1) = d(1)/2;
for i=2:N
    y(1,i) = (d(i) - md(1,i-1)*y(1,i-1))/p(1,i);
end
x = zeros(1,N);
x(1,N) = y(1,N);
for i=N-1:-1:1
    x(1,i) = y(1,i) - q(1,i)*x(1,i+1);
end
fprintf('M为三对角方程的解\n');
M=x;
fprintf('\n');
syms t;
digits(m);
for i=1:N-1
    pp(i)=M(i)*(X(i+1)-t)^3/(6*h(i))+M(i+1)*(t-X(i))^3/(6*h(i))+(Y(i)-M(i)*h(i)^2/6)*(X(i+1)-t)/h(i)+(Y(i+1)-M(i+1)*h(i)^2/6)*(t-X(i))/h(i);
    pp(i) = simplify(pp(i));
    coeff = sym2poly(pp(i));
    if length(coeff)~= 4
        tt=coeff(1:3);
        coeff(1:4)=0;
        coeff(2:4)=tt;
    end
    if x0>X(i) & x0<X(i+1)
        L=i;
        y0=coeff(1)*x0^3+coeff(2)*x0^2+coeff(3)*x0+coeff(4);
    end
    val = X(i):interval:X(i+1);
    for k=1:length(val)
        fval(k)=coeff(1)*val(k)^3+coeff(2)*val(k)^2+coeff(3)*val(k)+coeff(4);
    end
    if mod(i,2)==1
        plot(val,fval,'r+')
    else
        plot(val,fval,'b.')
    end
    hold on
    clear val fval
    ans = sym(coeff, 'd');
    ans = poly2sym(ans, 't');
    fprintf('在区间[%f,%f]内\n', X(i), X(i+1));
    fprintf('三次样条函数S(%d)=', i);
    pretty(ans);
end
fprintf('x0所在区间为[%f,%f]\n', X(L), X(L+1));
fprintf('函数在插值点x0=%f的值为\n', x0);
y0

end
