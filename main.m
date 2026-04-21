clear, clc;
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

m = 10;   % kg
v = 10;   % m/s
g = 9.81; % m/s^2
endingH = (g*yPeaks(1) - .5*(v)^2)/g;

finalKE = .5*m*(v)^2;
finalPE = m*g*yPeaks(1);

% Start by getting five points from user input, each representing a 
% peak/starting point. Set the final height to the highest peak minus four
% so that the final velocity is less than 10 m/s 
y0(1)=yPeaks(1);
y0(3)=yPeaks(2);
y0(5)=yPeaks(3);
y0(7)=yPeaks(4);
y0(9)=yPeaks(5);
y0(11)=endingH;

% Then, generate x coordinates. The peaks will be located every 60 meters
% (to make sure the coaster is at least 300 m). In between each peak, we
% have a random x-coordinate that is between 20 & 40 m from each peak.
x0=[0,randi([20,40]),60,randi([80,100]),120,randi([140,160]),180,randi([200,220]),240,randi([260,280]),300];

% Now, we generate troughs. The height of the troughs is located at a
% random location that is below the lower of the two peaks that it sits
% between.
y0(2)=randi([0,min([y0(1),y0(3)])]);
y0(4)=randi([0,min([y0(3),y0(5)])]);
y0(6)=randi([0,min([y0(5),y0(7)])]);
y0(8)=randi([0,min([y0(7),y0(9)])]);
y0(10)=randi([0,min([y0(9),y0(11)])]);

% Now, we use a assume that the rollercoaster path p(x) is modeled by a
% polynomial. As I will describe, we have 6 conditions that need to be met,
% so we will use a 5th order polynomial. For three conditions, the
% polynomial must equal the values of y at the values of x. For the other
% three, the derivative of y with respect to x at these three points must
% be zero. Through my calculations, I have determined that the polynomial 
% will take the form:
% y(x) = int_0^x((x-x0(1))(x-x0(2))(x-x0(3))(Ax+b))dx+y0(1)
%      = (A)/5.*x.^5+(b-A.*(x0(1)+x0(2)+x0(3))/4.*x.^4+(A*(x0(1)*x0(2)+x0(2)*x0(3)+x0(3)*x0(1))-b*(x0(1)+x0(2)+x0(3)))/3.*x.^3+(b*(x0(1)*x0(2)+x0(2)*x0(3)+x0(3)*x0(1))-A*x0(1)*x0(2)*x0(3))/2.*x.^2-b*x0(1)*x0(2)*x0(3).*x+y0(1)
% Also through my own calculations, I have found a pair of equations that
% can be used to solve for A and B: They are as follows
% Equation 1:
% ((A)/5.*x0(2).^5+(b-A.*(x0(1)+x0(2)+x0(3))/4.*x0(2).^4+(A*(x0(1)*x0(2)+x0(2)*x0(3)+x0(3)*x0(1))-b*(x0(1)+x0(2)+x0(3)))/3.*x0(2).^3+(b*(x0(1)*x0(2)+x0(2)*x0(3)+x0(3)*x0(1))-A*x0(1)*x0(2)*x0(3))/2.*x0(2).^2-b*x0(1)*x0(2)*x0(3).*x0(2))-((A)/5.*x0(1).^5+(b-A.*(x0(1)+x0(2)+x0(3))/4.*x0(1).^4+(A*(x0(1)*x0(2)+x0(2)*x0(3)+x0(3)*x0(1))-b)/3.*x0(1).^3+(b*(x0(1)*x0(2)+x0(2)*x0(3)+x0(3)*x0(1))-A*x0(1)*x0(2)*x0(3))/2.*x0(1).^2-b*x0(1)*x0(2)*x0(3).*x0(1))+y0(1)-y0(2)=0
% Equation 2:
% ((A)/5.*x0(3).^5+(b-A.*(x0(1)+x0(2)+x0(3))/4.*x0(3).^4+(A*(x0(1)*x0(2)+x0(2)*x0(3)+x0(3)*x0(1))-b*(x0(1)+x0(2)+x0(3)))/3.*x0(3).^3+(b*(x0(1)*x0(2)+x0(2)*x0(3)+x0(3)*x0(1))-A*x0(1)*x0(2)*x0(3))/2.*x0(3).^2-b*x0(1)*x0(2)*x0(3).*x0(3))-((A)/5.*x0(2).^5+(b-A.*(x0(1)+x0(2)+x0(3))/4.*x0(2).^4+(A*(x0(1)*x0(2)+x0(2)*x0(3)+x0(3)*x0(1))-b)/3.*x0(2).^3+(b*(x0(1)*x0(2)+x0(2)*x0(3)+x0(3)*x0(1))-A*x0(1)*x0(2)*x0(3))/2.*x0(2).^2-b*x0(1)*x0(2)*x0(3).*x0(2))+y0(2)-y0(3)=0
% This process must be carried out between each peak. We will make a for
% loop to do this.

for j=1:2:9

    syms A b;
    
    % Helpful intermediate sums to keep the equations readable
    S1 = x0(j) + x0(j+1) + x0(j+2);
    S2 = x0(j)*x0(j+1) + x0(j+1)*x0(j+2) + x0(j+2)*x0(j);
    S3 = x0(j)*x0(j+1)*x0(j+2);
    
    % Define the indefinite integral WITHOUT the constant of integration
    F = @(x, A_sym, b_sym) (A_sym/5).*x.^5 + (b_sym-A_sym*S1)/4.*x.^4 + (A_sym*S2-b_sym*S1)/3.*x.^3 + (b_sym*S2-A_sym*S3)/2.*x.^2 - b_sym*S3.*x;
    
    % The true expression shifts the curve so it evaluates exactly to y0(j) when x = x0(j)
    expr = @(x, A_sym, b_sym) F(x, A_sym, b_sym) - F(x0(j), A_sym, b_sym) + y0(j);
    
    % Set up the equations: y(x0(j+1)) = y0(j+1) & y(x0(j+2)) = y0(j+2)
    eq1 = expr(x0(j+1),A,b) == y0(j+1);
    eq2 = expr(x0(j+2),A,b) == y0(j+2);
    
    % Now, we need to solve the system of equations for A & b
    sol = solve([eq1, eq2], [A, b]);
    A_val = double(sol.A);
    b_val = double(sol.b);
    
    % Calculate the correct constant term (the y-intercept) for the polynomial array
    C_val = y0(j) - F(x0(j), A_val, b_val);
    
    % Create the polynomial vector using the clean S1, S2, S3 variables
    p = [A_val/5, (b_val-A_val*S1)/4, (A_val*S2-b_val*S1)/3, (b_val*S2-A_val*S3)/2, -b_val*S3, C_val];
    
    % Create space and evaluate
    x_space(j,:) = linspace(x0(j),x0(j+2),100);
    y(j,:)=polyval(p,x_space(j,:));

end
figure;
% Finally, we plot all of the polynomials as a piecewise function
plot(x0,y0,'*',x_space(1,:),y(1,:),x_space(3,:),y(3,:),x_space(5,:),y(5,:),x_space(7,:),y(7,:),x_space(9,:),y(9,:))


% plot(xPeaks, yPeaks, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Horizontal Distance (m)');
ylabel('Height (m)');
title('Roller Coaster Peak Inputs');
legend('Peak Heights', 'Location', 'best');
