function [A]=range0toN(A,range)
if nargin<2
    a=0;
    b=255;
else
    b=range(2);
    a=range(1);
end
A(A>b)=b;
A(A<a)=a;