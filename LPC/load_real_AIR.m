function [h_air_1, h_air_2] = load_real_AIR(fs)

airpar.fs = fs;
airpar.rir_type = 1; % Binaural AIR
airpar.room = 5; % Stairwell
airpar.head = 1; % HATS included
airpar.rir_no = 2; % 2m distance
airpar.azimuth = 90; % Front
airpar.channel = 1;
[h_air_1,air_info] = load_air(airpar);
h_air_1 = h_air_1';
airpar.channel = 0;
[h_air_2,air_info] = load_air(airpar);
h_air_2 = h_air_2';

end