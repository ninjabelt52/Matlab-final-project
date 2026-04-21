clc;
clear;
close all;

disp('Enter the heights of the five roller coaster peaks.');

yPeaks = zeros(1,5);

for i = 1:5
    yPeaks(i) = input(['Enter height for peak ', num2str(i), ': ']);
    
    while yPeaks(i) <= 0
        disp('Invalid input. Height must be greater than 0.');
        yPeaks(i) = input(['Enter height for peak ', num2str(i), ': ']);
    end
end

xPeaks(1) = randi([0,20]);

for i = 2:5
    xPeaks(i) = xPeaks(i-1) + randi([20,60]);
end

figure;
plot(xPeaks, yPeaks, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Horizontal Distance (m)');
ylabel('Height (m)');
title('Roller Coaster Peak Inputs');
legend('Peak Heights', 'Location', 'best');