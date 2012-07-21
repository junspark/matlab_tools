function [text_final, start_angle, finish_angle] = popla_convert(file_name)
%
%   USAGE:
%
%   [text_final, start_angle, finish_angle] = popla_convert(file_name)
%
%   INPUT:
%
%   file_name   is a column vector of strings,
%               it contains file names of GADDS data files
%
%   OUTPUT:
%
%   text_final  is a column vector,
%               it contains data obtained from GADDS in 5 degree increments
%               in azimuth and elevation
%
%   start_angle is a column vector,
%               it contains starting angular positiona (in degrees) of PF
%               coverage
%
%   finish_angle is a column vector,
%               it contains angular positions (in degrees) at which PF
%               coverage ends
%

fid = fopen(file_name);
white_space = [];
white_space_row = [];
text_final = [];

for i = 1:78
    text_line = fgetl(fid);
    if ~ischar(text_line)
        break
    end
    if i > 2
        white_space_row = isspace(text_line);
        white_space = [white_space; white_space_row(2:end)];
    end
end
fclose(fid);

text_digit = textread(file_name, '%1d', 'headerlines', 2);
white_space = reshape(white_space', prod(size(white_space)), 1);

j = 1;
for i = 1:length(white_space)
    if white_space(i) == 1
        text_col(i, 1) = 0;
    else
        text_col(i, 1) = text_digit(j);
        j = j + 1;
    end
end

j = 1;
for i = 1:4:length(text_col)
    text(j, 1) = 1000 * text_col(i) + 100 * text_col(i + 1) + 10 * text_col(i + 2) + text_col(i+3);
    j = j + 1;
end

format compact;
text = reshape(text, [18, 76])';
zerochk = find(sum(text, 2) == 0);

%fprintf('\nElements removed from GADDS file:\n\n');
for i = 1:size(text, 1)
    if ~isempty(zerochk)
        if i ~= zerochk
            text_final = [text_final; text(i, :)];
        else
            %% disp(text(i, :));
        end
    else
        text_final = [text_final; text(i, :)];
    end
end

text_final = reshape(text_final', [prod(size(text_final)), 1]);

if ~isempty(zerochk)
    i = 1;
    if zerochk(1) ~= 1
        text_final = text_final(72:end);
        start_angle = 0;
    else
        while i == zerochk(i)
            i = i + 1;
        end
        start_angle = (zerochk(i - 1)) * 1.25;
    end
    finish_angle = 90 - ((76 - zerochk(i) + 1) * 1.25);
else
    text_final = text_final(72:end);
    start_angle  = 0;
    finish_angle = 90;
end

text_final = text_final/100;