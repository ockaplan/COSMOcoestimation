function arcsec = rad2arcsec(rad)
%arcsec2rad - converts an angle in arcsecs to radians
%inputs: arcsec - angle in arcsecs
%outputs: theta - angle in rads
    arcsec = (rad .* 180 .* 3600) ./ pi;
end