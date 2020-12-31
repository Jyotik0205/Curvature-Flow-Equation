function message = Main
  assignin('base','Start',@Start);
  assignin('base','XDOT',@XDOT);
  assignin('base','PLOTN',@PLOTN);
  assignin('base','central_diff',@central_diff);
  message='Done importing functions to workspace';
end

function Start(M, N, t0, tf,c)
    % N  = #time steps
    % M  = #spatial discretization
    % t0 = initial time (t0 = 0)
    % tf = final time (tf = 4)
  % Curve 1
  if(c==1)
    x1_ = @(a) (4+cos(3*a)).*cos(a);
    x2_ = @(a) (4+cos(3*a)).*sin(a);
  end
   % Curve 2
  if(c==2)
     x1_ = @(a) (3+cos(5*a)).*cos(a);
     x2_ = @(a) (3+cos(5*a)).*sin(a);
  end
   % Curve 3
  if(c==3)
     x1_ = @(a) 3*cos(a);
     x2_ = @(a) 5*sin(a);
  end
    h = 2*pi/M;
    a = linspace(0, 2*pi-h, M);
    x = [x1_(a); x2_(a)];

    delx = (tf - t0)/N;

    for i = 0 : N
        x_dot = XDOT(x);
        if rem(i, N/(2*(tf-t0))) == 0
            PLOTN(x);
            pause(1);
        end  
        x = x + x_dot * delx;
    end
end

function x_dot = XDOT(x, M)
    % Derivates and double derivatives of x and y
    x_1a = central_diff(x(1,:));
    x_2a = central_diff(x(2,:));
    x_1aa = central_diff(x_1a);
    x_2aa = central_diff(x_2a);
    x_a = sqrt(x_1a.^2 + x_2a.^2);
    % Kappa calculation by formula
    kappa = (x_2a .* x_1aa - x_1a .* x_2aa) ./ x_a.^3;
    %Normals
    n1 = x_2a ./ x_a;
    n2 = -x_1a ./ x_a;
    %delx and dely
    x_dot1 = n1 .* kappa;
    x_dot2 = n2 .* kappa;
    x_dot = [x_dot1; x_dot2];
end

function [ result ] = central_diff( in )
% Computes the derivatives using central difference approximation
in=in';
% res = (in(2+in)-f(2-in))./(2*h);
M = size(in,1);
 h = (2*pi)/M;
in = [in(end);in;in(1)];
for i = 1: M
    result(i,1)=(in(i+2)-in(i))/2/h;
end

result=result';
end

function PLOTN(x)
    plot(x(1,:), x(2,:));
    axis equal
    hold on
   %hold off
end


 

