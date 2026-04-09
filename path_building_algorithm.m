clear;clc;

% Start by getting two points from user input, then generate a third point
% that has x between them and y below them. For now, the x-endpoints will
% be set at 0 & 60.
y0(1)=input('Please enter the starting height of the rollercoaster. Enter this as a scalar greater than 0.');
y0(3)=input('Please enter the height of the first peak of the rollercoaster. Enter this as a scalar greater than 0.');
x0=[0,randi([20,40]),60];
y0(2)=randi([0,min(y0)]);

% Now, we use a assume that the rollercoaster path p(x) is modeled by a
% polynomial. As I will describe, we have 6 conditions that need to be met,
% so we will use a 5th order polynomial. For three conditions, the
% polynomial must equal the values of y at the values of x. For the other
% three, the derivative of y with respect to x at these three points must
% be zero. Through my calculations, I have determined that the polynomial 
% will take the form:
% y(x) = int_0^x(x(x-x0(2))(x-60)(Ax+b))dx+y0(1)
% Also through my own calculations, I have found a pair of equations that
% can be used to solve for A and B: They are as follows
% Equation 1:
% A.*x0(2).^5/5+(b.*x0(2).^4-A.*x0(2).^5-60*A.*x0(2).^4)/4+(60*A.*x0(2).^4-b.*x0(2).^4-60*b.*x0(2).^3)/3+(60*b.*x0(2).^3)/2-y0(2)+y0(1)=0
% Equation 2:
% (A.*x0(3).^5/5+(b-A.*x0(2)-60*A).*x0(3).^4/4+(60*A.*x0(2)-b.*x0(2)-60*b).*x0(3)^3/3+(60*b.*x0(2)).*x0(3)^2/2)-(A.*x0(2).^5/5+(b.*x0(2).^4-A.*x0(2).^5-60*A.*x0(2).^4)/4+(60*A.*x0(2).^4-b.*x0(2).^4-60*b.*x0(2).^3)/3+(60*b.*x0(2).^3)/2)-y0(3)+y0(2)=0

syms A b;
eq1 = A*x0(2)^5/5 + (b*x0(2)^4 - A*x0(2)^5 - 60*A*x0(2)^4)/4 + (60*A*x0(2)^4 - b*x0(2)^4 - 60*b*x0(2)^3)/3 + (60*b*x0(2)^3)/2 - y0(2) + y0(1) == 0
eq2 = (A*x0(3)^5/5 + (b - A*x0(2) - 60*A)*x0(3)^4/4 + (60*A*x0(2) - b*x0(2) - 60*b)*x0(3)^3/3 + (60*b*x0(2))*x0(3)^2/2) - ...
      (A*x0(2)^5/5 + (b*x0(2)^4 - A*x0(2)^5 - 60*A*x0(2)^4)/4 + (60*A*x0(2)^4 - b*x0(2)^4 - 60*b*x0(2)^3)/3 + (60*b*x0(2)^3)/2) - y0(3) + y0(2) == 0

% Now, we need to solve the system of equations for A & b
sol = solve([eq1, eq2], [A, b]);
A_val = double(sol.A);
b_val = double(sol.b);

% Now, we create the polynomial and plot it alongside the user-requested
% points
p = [A_val/5,(b_val-A_val.*x0(2)-60.*A_val)/4,(60.*A_val.*x0(2)-b_val.*x0(2)-60.*b_val)/3,(60.*b_val.*x0(2))/2,0,y0(1)];
x_space = linspace(0,60,100);
y=polyval(p,x_space);
plot(x_space,y,x0,y0,'*')