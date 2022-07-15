% Abdelrahman Alaa
% Supervised by: Dr. Hesham Omran
% Integrated Circuits Laboratory (ICL)
% Electronics and Communications Eng. Dept.
% Faculty of Engineering
% Ain Shams University, Cairo, Egypt
%
% Version Date: 20200319
function GM = aaLookupGM(LUTV, L, VGS, VDS, VSB, TEMP)
% Assume LUT has a uniform spacing VGS Vector
step = LUTV.VGS(2) - LUTV.VGS(1);
% Get VGS values before and after VGS of query point (ends of interval)
VGS1 = LUTV.VGS(1) + floor((VGS - LUTV.VGS(1)) / step) * step;
VGS2 = LUTV.VGS(1) + ceil((VGS - LUTV.VGS(1)) / step) * step;
% Distance from VGS of query point to start of the interval
dVGS = VGS - VGS1;
% Function (Current) Value at start & end of the interval
ID1 = LUTV.ID(L, VGS1, VDS, VSB, TEMP);
ID2 = LUTV.ID(L, VGS2, VDS, VSB, TEMP);
% Slope (Transconductance) Value at start & end of the interval
GM1 = LUTV.GM(L, VGS1, VDS, VSB, TEMP);
GM2 = LUTV.GM(L, VGS2, VDS, VSB, TEMP);
GM = (6 * step * dVGS - 6 * dVGS.^2) ./ step^3 .* ID2 ...
    + (-6 * step * dVGS + 6 * dVGS.^2) ./ step^3 .* ID1 ...
    + (3 * dVGS.^2 - 2 * step * dVGS) ./ step^2 .* GM2 ...
    + (3 * dVGS.^2 - 4 * step * dVGS + step^2) ./ step^2 .* GM1;
