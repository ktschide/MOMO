function [row, col] = UTMlatlon2pix(R, lat, lon)

dx = R(2,1);
dy = R(1,2);

xll = R(3,1);
yll = R(3,2);

row = (lat-yll)/dy;
col = (lon-xll)/dx;
end