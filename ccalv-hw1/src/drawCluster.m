function cellCenters = drawCluster( centerPosition, iValue, jValue, cellRadius, clusterNumber)
% Cameron Calv


N = iValue^2 + iValue*jValue + jValue^2;
%Do the I and J values correspond to a proper clustering?
cellLetters = [];
if (N == 3)
    cellLetters = ['A', 'B', 'C'];
elseif (N == 4)
    cellLetters = ['A', 'B', 'C', 'D'];
elseif (N == 7)
    cellLetters = ['A', 'B', 'C', 'D', 'E', 'F', 'G'];
else
    return;
end

incircleRadius = (sqrt(3)/2)*cellRadius;
cellCenters = [];

for n=1:N
    if n>1
        nextCenter = centerPosition + ...
            2*incircleRadius*(cos((n-1)*(pi/3)-pi/6)...
            +j*sin((n-1)*(pi/3)-pi/6));
        drawCell( nextCenter, cellRadius, cellLetters(n)+"_"+num2str(clusterNumber) );
        cellCenters = [cellCenters nextCenter];
    else
        drawCell( centerPosition, cellRadius, cellLetters(n)+"_"+num2str(clusterNumber) );
        cellCenters = [cellCenters centerPosition];
    end
end
end
