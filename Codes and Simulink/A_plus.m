% System parameters
clc;clear,clf;
m = 0.03; % unbalanced mass
e = 0.02; % eccentricity of unbalanced mass
w =2500*2*pi/60 ; % angular velocity of motor (rad/sec)
F_amp = m*e*w^2 ; % Amplitude of the excitation force (N)
m1 = 0.55; % Mass of the pump (kg)
m2 = 0.125*2;  % Mass of the absorber (kg)
k1 = 38000; % Stiffness of the pump (N/m)
k2 = 17000;  % Stiffness of the absorber (N/m)
w_n = (k1/m1)^0.5;  % natural frequency of the system
w_a = (k2/m2)^0.5;  % natural frequency of the absorber
zeta=0.05;   % damping factor
c_1=zeta*2*m1*w_n;  % damping coeffecient

a=(k1+k2)/m1;b=-k2/m1;c=-k2/m2;d=k2/m2;e=a+d;f=a*d-b*c;

w_n1=((e-(e^2-4*f)^0.5)/2)^0.5;
w_n2=((e+(e^2-4*f)^0.5)/2)^0.5;

a1=(a-w_n1^2)/-b;
b1=(a-w_n2^2)/-b;

x_1=((k2-w^2*m2)*F_amp)/((k1+k2-w^2*m1)*(k2-w^2*m2)-k2^2);
x_2=(k2*F_amp)/((k1+k2-w^2*m1)*(k2-w^2*m2)-k2^2);
B2=((x_2-a1*x_1)*w)/((b1-a1)*w_n2);
B1=-(x_1*w+B2*w_n2)/w_n1;
r=w/(k1/m1)^0.5;
x_st=F_amp/k1;
x_o=(x_st)/((1-r^2)^2+(2*r*zeta)^2)^0.5;
if zeta==0
    t=0:0.0001:25;
    if r==0
        x=x_st/2*w_n*t;
    else
        x=x_o*(sin(w*t)-r*sin(w_n*t));
    end
    x1=x_1*(sin(w*t))+B1*sin(w_n1*t)+B2*sin(w_n2*t);
    x2=x_2*(sin(w*t))+a1*B1*sin(w_n1*t)+b1*B2*sin(w_n2*t);
    figure(1)
    plot(t,x)
    grid on
    xlabel('time (sec)');
    ylabel('amplitude (m)');
    title('Machine Response without absorber' );
    
    figure(2)
    plot(t,x1)
    grid on
    xlabel('time (sec)');
    ylabel('amplitude (m)');
    title('Machine Response with absorber' );
    
    figure(3)
    plot(t,x2)
    grid on
    xlabel('time (sec)');
    ylabel('amplitude (m)');
    title('Response of absorber' );
elseif zeta<0
    t=0:0.0001:0.2;
    fprintf("wrong value for zeta")
elseif zeta~=0
    t=0:0.0001:0.2;
    x=x_o*(sin(w*t));
    x1=x_1*(sin(w*t));  
    x2=x_2*(sin(w*t));

    figure(1)
    plot(t,x)
    grid on
    xlabel('time (sec)');
    ylabel('amplitude (m)');
    title('Machine Response without absorber' );

    figure(2)
    plot(t,x1)
    grid on
    xlabel('time (sec)');
    ylabel('amplitude (m)');
    title('Machine Response with absorber' );
    
    figure(3)
    plot(t,x2)
    grid on
    xlabel('time (sec)');
    ylabel('amplitude (m)');
    title('Response of absorber' );

end



% Constants

ws = sqrt(k1 / m1); % Natural frequency of the main system (rad/s)

% Mass ratios to analyze
w_wa = linspace(0, 2, 100000); % Frequency ratio (w/wa)

% Set up the figure


% Parameters for each mass ratio
m2 = 0.125*2;          % Mass of absorber
k2 = 17000;                  % Stiffness of the absorber (N/m)
wa = sqrt(k2 / m2);       % Absorber natural frequency (rad/s)

% Machine response: X1/Xst
X1_Xst = abs((1 - (w_wa).^2) ./ ...
             ((1 + k2/k1 - (w_wa * ws / wa).^2) .* ...
              (1 - (w_wa).^2) - k2/k1));

% Absorber response: X2/Xst
X2_Xst = abs(1 ./ ...
             ((1 + k2/k1 - (w_wa * ws / wa).^2) .* ...
              (1 - (w_wa).^2) - k2/k1));

% Plot machine response
figure(4)

subplot(1, 2,1);
plot(w_wa, X1_Xst, 'r', 'LineWidth', 1.5);
xlabel('Frequency ratio (\omega / \omega_a)');
ylabel('abs(X_1 / X_{st})');
title('Machine Response' );
grid on;
axis([0 2 0 20]);

% Plot absorber response
subplot(1, 2,2);
plot(w_wa, X2_Xst, 'b', 'LineWidth', 1.5);
xlabel('Frequency ratio (\omega / \omega_a)');
ylabel('abs(X_2 / X_{st})');
title('Absorber Response' );
grid on;
axis([0 2 0 20]);

Reg_1=abs((w-w_n1)/w*100);
Reg_2=abs((w-w_n2)/w*100);
reg=[Reg_1  Reg_2];

fprintf("\nThe natural frequencies is at least : \n +- %f %% from the operating speed\n",min(reg))

