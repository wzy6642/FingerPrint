# FingerPrint
指纹识别算法总结
spline3.m是三次样条插值函数，可以输出插值多项式，使用方法如下：

X = 0.001:3:6.001;

Y = sin(X)./X;

dY=[-0.000333333 0.16778];

m=5;

x0=1.5;

spline3(X,Y,dY,x0,m)

其中dY是通过mathematica求解得到，spline3.nb是对上述算法的验证文件
