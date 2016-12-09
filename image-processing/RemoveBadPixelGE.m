function im_out = RemoveBadPixelGE(im_in, b2, b1)
neg_pixel               = im_in < 0;

im_out      = im_in;
im_out(b2)  = (im_out(b2 - 2048) + im_out(b2 + 2048) + im_out(b2 - 1) + im_out(b2 + 1))/4;
im_out(neg_pixel | b1)  = 0;


%     if (BadFileContents[ElementNr]==2 && ElementNr > 2048 && ElementNr < (2047*2048)){
%                         Data[ElementNr] = (Data[ElementNr-1] + Data[ElementNr+1] + Data[ElementNr + 2048] + Data[ElementNr - 2048])/(4.0);
%                     }
%                     if (Data[ElementNr] < 0 || (BadFileContents[ElementNr])%2 == 1) {
%                         Data[ElementNr] = 0;
%                     }