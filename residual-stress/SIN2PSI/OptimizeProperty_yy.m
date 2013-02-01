function ydata = OptimizeProperty_yy(x, xdata)
E   = x(1);
nu  = x(2);

S1  = E/(1 + nu);
S2  = E/nu;

m   = xdata(:,1);
b   = xdata(:,2);

%%% FOR SIGXX
sigxx   = S1*m;

%%% FOR SIGYY
sigyy   = - S2*b - sigxx;

% ydata   = [sigxx; sigyy];
% ydata   = sigxx;
ydata   = sigyy;