function f=dcum(a,b,t)

rd=pi/180;
if (a>=1)
    f = 1-cos(rd*t);
else
    eps=1e-8;
    delx=1;
    x=2*rd*t;
    p=x;
	while (delx >= eps)
        y = a*sin(x)+.5*b*sin(2.*x);
        dx=.5*(y-x+p);
        x=x+dx;
        delx=abs(dx);
    end
	f = (2.*y+p)/pi;
end
