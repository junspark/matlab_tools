function CorrectedImageData = CorrectBadPixels(ImageData, BadPixelData)
% CorrectedImageData - output corrected image based on baspixel data for GE
% images
%
%   INPUT:
%
%   ImageData
%       raw image data
%
%   BadPixelData
%       list of bad pixel number
%
%   OUTPUT:
%
%   CorrectedImageData
%       Corrected image
%
%   For a given bad pixel, the four neighboring pixels are averged and the
%   average is assigned. If the bad pixel is an edge pixel, only the pixels
%   inside the image are avearged.


CorrectedImageData  = ImageData;
nbr1    = BadPixelData - 2048 - 1;
nbr2	= BadPixelData - 2048 + 1;
nbr3    = BadPixelData + 2048 - 1;
nbr4    = BadPixelData + 2048 + 1;

n   = length(BadPixelData);
for i = 1:1:n
    ct  = 0; p1  = 0; p2  = 0; p3  = 0; p4  = 0;
    if nbr1(i) <= 2048*2048 && nbr1(i) >= 1
        p1  = ImageData(nbr1(i));
        ct  = ct + 1;
    end
    if nbr2(i) <= 2048*2048 && nbr2(i) >= 1
        p2  = ImageData(nbr2(i));
        ct  = ct + 1;
    end
    if nbr3(i) <= 2048*2048 && nbr3(i) >= 1
        p3  = ImageData(nbr3(i));
        ct  = ct + 1;
    end
    if nbr4(i) <= 2048*2048 && nbr4(i) >= 1
        p4  = ImageData(nbr4(i));
        ct  = ct + 1;
    end
    p   = (p1 + p2 + p3 + p4)/ct;
    
    CorrectedImageData(BadPixelData(i)) = p;
end