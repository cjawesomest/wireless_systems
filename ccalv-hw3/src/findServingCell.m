function [cellNumber, tierNumber, center] = findServingCell( mobileLocation, cellCenters)

distance = inf;
center = [];
element_number = 0;
for c=1:numel(cellCenters)
    nextDist_real = real(mobileLocation)-real(cellCenters(c));
    nextDist_imag = imag(mobileLocation)-imag(cellCenters(c));
    nextDistance = sqrt(nextDist_real^2 + nextDist_imag^2);
    if(abs(nextDistance) <= abs(distance))
        distance = nextDistance;
        element_number = c;
        center = cellCenters(element_number);
    end
end
tierNumber = mod(element_number-1, size(cellCenters, 2))+1;
cellNumber = floor((element_number-1)/size(cellCenters, 2))+1;

end