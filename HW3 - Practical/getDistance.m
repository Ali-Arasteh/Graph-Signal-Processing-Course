function D = getDistance(lat1, lon1, lat2, lon2)
rad = pi / 180;
radius = 6371; %earth radius in kilometers
D = abs(acos(sin(lat2 * rad) * sin(lat1 * rad) + cos(lat2 * rad) * cos(lat1 * rad) * cos(lon2 * rad - lon1 * rad)) * radius); %result in Kilometers
end