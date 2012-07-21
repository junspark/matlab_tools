function [V, D] = makeDiscHarm(h1sip, l2ip, numh); 
% [V, D] = makeDiscHarm(h1sip, l2ip, numh);
  
  [V, D] = eigs(h1sip, l2ip, numh, 'sm');
  [D, ind] = sort(diag(D));
  V = V(:, ind);
  
  