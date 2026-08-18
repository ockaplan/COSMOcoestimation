function theta = arcsec2rad(arcsec)
%arcsec2rad - converts an angle in arcsecs to radians
%inputs: arcsec - angle in arcsecs
%outputs: theta - angle in rads
    theta = (arcsec * pi) / (180*3600);
end