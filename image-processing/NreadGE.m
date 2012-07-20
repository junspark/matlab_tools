function img = NreadGE(filename, frameno)
buffer  = 8192;

fp      = fopen(filename,'r','n');
offset  = buffer + (frameno-1)*2048*2048*2;

fseek(fp,offset,'bof');

img     = fread(fp,[2048 2048],'uint16');

fclose(fp);
