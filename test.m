clear all; close all; clc;
f = @(x) 1./(1+25*x.^2);
p0 = @(x) 1;
p1 = @(x) x;

alpha = sqrt(1/3);
g2pt = @(a, b, n, f) (b-a)/n/2*sum(f(((b-a)/n/2:(b-a)/n:b)-alpha*(b-a)/n/2) + f(((b-a)/n/2:(b-a)/n:b)+alpha*(b-a)/n/2));
g = @(x) f(x).*p1(x);

a = 0;
b = 1;
n = 128;
b1 = g2pt(a,b,n,g);
b0 = g2pt(a,b,n,f);
h = @(x) p0(x).*p0(x);
B11 = 1;
h = @(x) p0(x).*p1(x);
B12 = g2pt(a,b,n,h);
B21 = g2pt(a,b,n,h);
h = @(x) p1(x).*p1(x);
B22 = g2pt(a,b,n,h);
B = [B11, B12; B21, B22];
b = [b0;b1];
c = linsolve(B,b);
x = [a:0.01:1];
v = c(1) + c(2)*x;
plot(x,v,'r-', x,f(x),'o')