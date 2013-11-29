function [y_out,x_out,c_val]=subpixel_detection(c,s)

%first I check where the max is, and construct a set of numbers around the
%maximum.
%then I will use the spline function to contstruc the optimal set. 
%then I will use the fnval function to evaluate the function and look for
%the maximum using fminserach on the negative function

[max_c, imax] = max(c(:));
[ypeak, xpeak] = ind2sub(size(c),imax(1));

c_cut=c(ypeak-s:ypeak+s,xpeak-s:xpeak+s);
[j,i]=ind2sub([2*s+1,2*s+1],[1:(2*s+1)^2]);
xy=[j;i];
for n=1:length(xy)
    val(n)=c_cut(xy(1,n),xy(2,n));
end

st = tpaps(xy,val);
%fnplt(st);
xf = fminsearch(@(x) fun(x,st),[s+1;s+1]);
y_out=xf(1)+ypeak-s-1;
x_out=xf(2)+xpeak-s-1;
c_val=fnval(st,xf);
function y=fun(x,st)
y=-fnval(st,x);







