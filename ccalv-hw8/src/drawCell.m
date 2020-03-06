function drawCell( center, radius, cellName )
% drawCell.m
%     Kapil R. Dandekar
%     ECE-T512 - Wireless Systems Matlab Courseware
%        
%     This program draws a hexagonal cell with
% a label can be drawn in a matlab figure.
%
%  Note that the coordinate system is x+j*y
%    i.e. X is real axis, and Y is imaginary axis
%  Remember to be consistent with your units
%
%  Parameters:
%     center - Location of cell center specified as a complex
%         number: x_center + j*y_center
%     radius - Cell radius specified as a real number
%     cellName - String labeling the name of the cell - note that
%         LaTex formatting can be used in matlab to make text
%         look better

plot( center  + radius * exp(j*pi*(0:2:12)/6), 'k');
text( real(center), imag(center), cellName );