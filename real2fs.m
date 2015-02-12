function xyfs = real2fs(gridpt, rad)
tic
out1 = [(gridpt(1)*cos(rad)+gridpt(2)*sin(rad)) (gridpt(2)*cos(rad)-gridpt(1)*sin(rad))]
out1time = toc
tic
out2 = [dot(gridpt,[cos(rad) sin(rad)]), dot(gridpt,[cos(rad-pi/2) sin(rad-pi/2)])]
out2time = toc

xyfs = 0;

end