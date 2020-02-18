function cellCenters = drawManyClusters( centerPosition, iValue, jValue, cellRadius, numberOfClusters)
%Cameron Calv

cellCenters = [];
for c=1:numberOfClusters
    if (c == 1)
        clusterCenters = drawCluster( centerPosition, iValue, jValue, cellRadius, c);
    else
        direction = pi/6 + (c-1)*(pi/3);
        nextPosition = centerPosition + 2*(sqrt(3)/2)*cellRadius*iValue*(cos(direction)+j*sin(direction));
        direction = direction + pi/3;
        nextCenter = nextPosition + 2*(sqrt(3)/2)*cellRadius*jValue*(cos(direction)+j*sin(direction));
        clusterCenters = drawCluster(nextCenter, iValue, jValue, cellRadius, c);
    end
    cellCenters = [cellCenters; clusterCenters];
end

end