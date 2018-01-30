% Give an interger n and number of digit m
% returns m-digit string n padded with zero
%
% Updated: 2014.02.26 by Andrew Chuang, APS

function string = padzero(n,numdigitout);

if nargin == 0
    fprintf('\nUsage: str = padzero(integer_n,number_of_digit_output)\n\n');
    return;
end

if numdigitout<0 || rem(numdigitout,1)~=0
    fprintf('output number of digit must be a positve integer\n');
    return;
end

if n<0 || rem(n,1)~=0;
    fprintf('\nn must be a positive integer\n\n');
    return;
end

if n == 0
    numdigit = 1;
else
    numdigit = floor(log10(n)) + 1;  % this works for non-integer 
end

numdigitout = max(numdigitout,numdigit);
string=repmat('0',1,numdigitout);
string(end-numdigit+1:end)=num2str(n);
