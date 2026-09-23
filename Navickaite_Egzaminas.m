clc
clear all
close all

% Iejimai: x(n) ir x(n-1). 
% Du iejimai, 2 paslepti sluoksniai, 3 isejimai.
% 1 sl. tiesinis, 2 sl. tanh, 3 sl. sigmoide (pilnai sujungtas tinklas).
% Tikslo funkcija: xi = 0.5*e^2, isvestine dxi/de = e.
%cia formules DI sugeneravo kad butu kazkoks rezultatas matomas:
x = 0.1:1/22:1;
d1 = ((1 + 0.6*sin(2*pi*x/0.7)) + 0.3*sin(2*pi*x))/2;
d2 = (1 + sin(2*pi*x))/2;
d3 = (1 + cos(2*pi*x))/2;

% Pradines koeficientu reiksmes
w11_1 = randn*0.5; w12_1 = randn*0.5; b1_1 = randn*0.5;
w21_1 = randn*0.5; w22_1 = randn*0.5; b2_1 = randn*0.5;
w31_1 = randn*0.5; w32_1 = randn*0.5; b3_1 = randn*0.5;

w11_2 = randn*0.5; w12_2 = randn*0.5; w13_2 = randn*0.5; b1_2 = randn*0.5;
w21_2 = randn*0.5; w22_2 = randn*0.5; w23_2 = randn*0.5; b2_2 = randn*0.5;

w11_3 = randn*0.5; w12_3 = randn*0.5; b1_3 = randn*0.5;
w21_3 = randn*0.5; w22_3 = randn*0.5; b2_3 = randn*0.5;
w31_3 = randn*0.5; w32_3 = randn*0.5; b3_3 = randn*0.5;

LR = 0.1;
epochos = 5000;
kastai = zeros(1, epochos);

for iter = 1:epochos
    kastai_ep = 0;

    for i = 2:length(x) % pradedam nuo 2, nes reikia x(n-1)
        xn  = x(i);
        xn1 = x(i-1);

        % 1 sluoksnio pasvertos sumos
        v1_1 = w11_1*xn + w12_1*xn1 + b1_1;
        v2_1 = w21_1*xn + w22_1*xn1 + b2_1;
        v3_1 = w31_1*xn + w32_1*xn1 + b3_1;
        % 1 sluoksnio aktyvacija (tiesine)
        y1_1 = v1_1;
        y2_1 = v2_1;
        y3_1 = v3_1;

        % 2 sluoksnio pasvertos sumos
        v1_2 = w11_2*y1_1 + w12_2*y2_1 + w13_2*y3_1 + b1_2;
        v2_2 = w21_2*y1_1 + w22_2*y2_1 + w23_2*y3_1 + b2_2;
        % 2 sluoksnio aktyvacija (tanh)
        y1_2 = tanh(v1_2);
        y2_2 = tanh(v2_2);

        % 3 sluoksnio (isejimo) pasvertos sumos
        v1_3 = w11_3*y1_2 + w12_3*y2_2 + b1_3;
        v2_3 = w21_3*y1_2 + w22_3*y2_2 + b2_3;
        v3_3 = w31_3*y1_2 + w32_3*y2_2 + b3_3;
        % Tinklo atsakas (sigmoide)
        y1 = 1/(1+exp(-v1_3));
        y2 = 1/(1+exp(-v2_3));
        y3 = 1/(1+exp(-v3_3));

        % Klaidos ir tikslo funkcija xi = 0.5*(e1^2+e2^2+e3^2)
        e1 = d1(i) - y1;
        e2 = d2(i) - y2;
        e3 = d3(i) - y3;
        kastai_ep = kastai_ep + 0.5*(e1^2+e2^2+e3^2);

        % 3 (isejimo) sluoksnio gradientai (sigmoides isvestine y(1-y))
        delta1_3 = y1*(1-y1)*e1;
        delta2_3 = y2*(1-y2)*e2;
        delta3_3 = y3*(1-y3)*e3;
        % 2 sluoksnio gradientai (tanh isvestine 1-y^2)
        delta1_2 = (1-y1_2^2) * (delta1_3*w11_3 + delta2_3*w21_3 + delta3_3*w31_3);
        delta2_2 = (1-y2_2^2) * (delta1_3*w12_3 + delta2_3*w22_3 + delta3_3*w32_3);
        % 1 sluoksnio gradientai (tiesines isvestine = 1)
        delta1_1 = delta1_2*w11_2 + delta2_2*w21_2;
        delta2_1 = delta1_2*w12_2 + delta2_2*w22_2;
        delta3_1 = delta1_2*w13_2 + delta2_2*w23_2;

        % 3 (isejimo) sluoksnio koeficientu atnaujinimas
        w11_3 = w11_3 + LR*delta1_3*y1_2;
        w12_3 = w12_3 + LR*delta1_3*y2_2;
        b1_3  = b1_3  + LR*delta1_3;
        w21_3 = w21_3 + LR*delta2_3*y1_2;
        w22_3 = w22_3 + LR*delta2_3*y2_2;
        b2_3  = b2_3  + LR*delta2_3;
        w31_3 = w31_3 + LR*delta3_3*y1_2;
        w32_3 = w32_3 + LR*delta3_3*y2_2;
        b3_3  = b3_3  + LR*delta3_3;
        % 2 sluoksnio koeficientu atnaujinimas
        w11_2 = w11_2 + LR*delta1_2*y1_1;
        w12_2 = w12_2 + LR*delta1_2*y2_1;
        w13_2 = w13_2 + LR*delta1_2*y3_1;
        b1_2  = b1_2  + LR*delta1_2;
        w21_2 = w21_2 + LR*delta2_2*y1_1;
        w22_2 = w22_2 + LR*delta2_2*y2_1;
        w23_2 = w23_2 + LR*delta2_2*y3_1;
        b2_2  = b2_2  + LR*delta2_2;
        % 1 sluoksnio koeficientu atnaujinimas
        w11_1 = w11_1 + LR*delta1_1*xn;
        w12_1 = w12_1 + LR*delta1_1*xn1;
        b1_1  = b1_1  + LR*delta1_1;
        w21_1 = w21_1 + LR*delta2_1*xn;
        w22_1 = w22_1 + LR*delta2_1*xn1;
        b2_1  = b2_1  + LR*delta2_1;
        w31_1 = w31_1 + LR*delta3_1*xn;
        w32_1 = w32_1 + LR*delta3_1*xn1;
        b3_1  = b3_1  + LR*delta3_1;
    end
    kastai(iter) = kastai_ep/(length(x)-1);
end

% Atsakas mokymo duomenims, be koeficientu atnaujinimo
Y1 = zeros(size(d1));
Y2 = zeros(size(d2));
Y3 = zeros(size(d3));
for i = 2:length(x)
    xn  = x(i);
    xn1 = x(i-1);
    y1_1 = w11_1*xn + w12_1*xn1 + b1_1;
    y2_1 = w21_1*xn + w22_1*xn1 + b2_1;
    y3_1 = w31_1*xn + w32_1*xn1 + b3_1;
    y1_2 = tanh(w11_2*y1_1 + w12_2*y2_1 + w13_2*y3_1 + b1_2);
    y2_2 = tanh(w21_2*y1_1 + w22_2*y2_1 + w23_2*y3_1 + b2_2);
    Y1(i) = 1/(1+exp(-(w11_3*y1_2 + w12_3*y2_2 + b1_3)));
    Y2(i) = 1/(1+exp(-(w21_3*y1_2 + w22_3*y2_2 + b2_3)));
    Y3(i) = 1/(1+exp(-(w31_3*y1_2 + w32_3*y2_2 + b3_3)));
end

% Testavimas. 
x_test  = 0:1/44:1;
d1_test = ((1 + 0.6*sin(2*pi*x_test/0.7)) + 0.3*sin(2*pi*x_test))/2;
d2_test = (1 + sin(2*pi*x_test))/2;
d3_test = (1 + cos(2*pi*x_test))/2;
Y1_test = zeros(size(x_test));
Y2_test = zeros(size(x_test));
Y3_test = zeros(size(x_test));
for i = 2:length(x_test)
    xn  = x_test(i);
    xn1 = x_test(i-1);
    y1_1 = w11_1*xn + w12_1*xn1 + b1_1;
    y2_1 = w21_1*xn + w22_1*xn1 + b2_1;
    y3_1 = w31_1*xn + w32_1*xn1 + b3_1;
    y1_2 = tanh(w11_2*y1_1 + w12_2*y2_1 + w13_2*y3_1 + b1_2);
    y2_2 = tanh(w21_2*y1_1 + w22_2*y2_1 + w23_2*y3_1 + b2_2);
    Y1_test(i) = 1/(1+exp(-(w11_3*y1_2 + w12_3*y2_2 + b1_3)));
    Y2_test(i) = 1/(1+exp(-(w21_3*y1_2 + w22_3*y2_2 + b2_3)));
    Y3_test(i) = 1/(1+exp(-(w31_3*y1_2 + w32_3*y2_2 + b3_3)));
end

fprintf('mokymo |e1| = %.4f, |e2| = %.4f, |e3| = %.4f\n', ...
    mean(abs(d1(2:end)-Y1(2:end))), mean(abs(d2(2:end)-Y2(2:end))), mean(abs(d3(2:end)-Y3(2:end))));
fprintf('testo  |e1| = %.4f, |e2| = %.4f, |e3| = %.4f\n', ...
    mean(abs(d1_test(2:end)-Y1_test(2:end))), mean(abs(d2_test(2:end)-Y2_test(2:end))), mean(abs(d3_test(2:end)-Y3_test(2:end))));

figure
subplot(1,3,1)
plot(d1(2:end), Y1(2:end), 'b.'); hold on
plot(d1_test(2:end), Y1_test(2:end), 'r.'); hold off
xlabel('norimas d1'); ylabel('tinklo y1');
legend('mokymas','testavimas'); title('1 isejimas'); grid on

subplot(1,3,2)
plot(d2(2:end), Y2(2:end), 'b.'); hold on
plot(d2_test(2:end), Y2_test(2:end), 'r.'); hold off
xlabel('norimas d2'); ylabel('tinklo y2');
legend('mokymas','testavimas'); title('2 isejimas'); grid on

subplot(1,3,3)
plot(d3(2:end), Y3(2:end), 'b.'); hold on
plot(d3_test(2:end), Y3_test(2:end), 'r.'); hold off
xlabel('norimas d3'); ylabel('tinklo y3');
legend('mokymas','testavimas'); title('3 isejimas'); grid on

% Papildomai - atsakas prieš x, kaip anksciau darei
figure
plot(x,d1,'b*'); hold on; plot(x,Y1,'r'); hold off
legend('Norimas atsakas (d1)','Tinklo atsakas (Y1)')

figure
plot(x,d2,'b*'); hold on; plot(x,Y2,'r'); hold off
legend('Norimas atsakas (d2)','Tinklo atsakas (Y2)')

figure
plot(x,d3,'b*'); hold on; plot(x,Y3,'r'); hold off
legend('Norimas atsakas (d3)','Tinklo atsakas (Y3)')
