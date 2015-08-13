function GM = ReadGrainMap( name ) 

fid = fopen(name);
tline = fgetl(fid);
counter = 1;

while ischar(tline)
    
    A = sscanf( tline , '%f', 10)';
    GM{1,counter} = A(1);
    GM{2,counter} = A(2:4);
    GM{3,counter} = A(5:7);
    GM{4,counter} = A(8);
    GM{5,counter} = A(9);
    GM{6,counter} = A(10);
    B = sscanf(tline,'%f',10+2*A(10) );
    GM{7,counter} = B(11:11+A(10)-1);
    GM{8,counter} = B(11+A(10):end);
    tline = fgetl(fid);
    counter = counter+1;
    
end

fclose(fid);