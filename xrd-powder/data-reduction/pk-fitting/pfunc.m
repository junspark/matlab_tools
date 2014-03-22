function y = pfunc(p,x)
p   = p(:);
x   = x(:);
y   = x.*0;

numpk   = (length(p) - 2)/4;
for i = 1:1:numpk
    ji  = 4*(i - 1) + 1;
    jf  = 4*i;
    ppk = p(ji:jf);
    ypk = pkpseudoVoigt(ppk,x);
    y   = y + ypk;
end
pbkg    = p((numpk*4+1):end);
ybkg    = polyval(pbkg,x);

y   = y + ybkg;