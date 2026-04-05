clear, clc;

startHeight = input("Please enter your starting height");

m = 10;   % kg
v = 10;   % m/s
g = 9.81; % m/s^2
endingH = (g*startHeight - .5*(v)^2)/g;

finalKE = .5*m*(v)^2;
finalPE = m*g*startHeight;