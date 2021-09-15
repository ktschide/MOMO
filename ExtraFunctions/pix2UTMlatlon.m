function [lat, lon] = pix2UTMlatlon(R, row, col)

dx = R(2,1);
dy = R(1,2);

xll = R(3,1);
yll = R(3,2);

lat = row*dy+yll;
lon = col*dx+xll;
end