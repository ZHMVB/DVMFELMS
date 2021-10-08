%% 2020-7-21 NSP_ANC

clear;
%% generate chaotic nois

M = 40000;
T = 1;

x = sin(2*pi*500*(1:M)/8000);
x = x';

%% generate d(n)
N = 10; 
% Pz主路径
Pz = zeros(N,1);
Pz(6) = 1; Pz(7) = 0.8; Pz(8) = 0.3;Pz(9)=0.4;

% Sz 次级路径
N1 = 3;
N2 = 3;
S1 = zeros(N1,1);
S2 = zeros(N2,1);
S1(1) = 1; S1(2) = -0.6; S1(3) = 0.05;
S2(1) = 1; S2(2) = 0.2; S2(3) = 0.05;

d = zeros(M,1);
x_input = zeros(N,1);
for i = 1:M
    x_input = [x(i);x_input(1:end-1)];
    d(i) = x_input(6) + 0.8 * x_input(7) + 0.3 * x_input(8) + 0.4 * x_input(9) + 0.2 * x_input(6) * x_input(7) - 0.3 * x_input(6) * x_input(8) + 0.4 * x_input(6) * x_input(9);
end

SNR = 30; % dB
SNR_B = 10^(SNR/20);
aaa = randn(M,1);
e_index = 1:floor(M/5);
EEE1 = sum(aaa(e_index).^2);
EEE2 = sum(d(e_index).^2);
SNR_B = sqrt(EEE2/ EEE1) / SNR_B ;

dd = d + SNR_B * randn(M,1);

EEEd = dd(1:floor(M/5));
EEEd = EEEd.^2;
EEEd = sum(EEEd(1:floor(M/5)))/floor(M/5);
EEE = 10*log10(EEEd);

dd = d;

%% FELMS 算法

miu1 = 1e-1;
S1_E = flipud(S1);
S2_E = flipud(S2);

E_main = zeros(M,1);
Beita = N1+N2-2;

for k = 1:T
    
    x_index = 1:M;
    xx = x(x_index);
    d = [dd(x_index(1:M/2));-dd(x_index(M/2+1:M))];
    d = awgn(d,SNR);
    
    x_input = zeros(N+N,1);
    
    X1n = zeros(1,1);
    X1 = zeros(N,1);
    X2 = zeros(N,1);
    W1 = zeros(N,1);
   
    Yn = 0;
    Y = zeros(N1+N2,1);    
    FE = zeros(N1+N2,1);
    M1 = zeros(N1+N2,1);
    M2 = zeros(N1+N2,1);

    i = 1;
    e = zeros(M,1);
    e_main = zeros(M,1);
    ee = EEEd ;
    ddd = EEEd ; 

    for i = 1:M
        
        %产生x_input 信号
        x_input = [x(i);x_input(1:end-1)];
        
        %产生Volterra X信号
        X1 = x_input(1:N);
        X2 = x_input(Beita+1:Beita+N);
        
        %产生输出信号
        Yn = X1' * W1;
        Y = [Yn;Y(1:end-1)];
        [Y_s,M1] = f_Sz(Y,M1,S1,S2,N1,N2);         
        e(i) = d(i) + Y_s;        
        
        ee = 0.95 * ee + 0.05 * e(i)^2;
        ddd = 0.95 * ddd + 0.05 * d(i)^2;
        e_main(i) = ee/ddd;
        
        % 滤波误差信号
        FE = [e(i);FE(1:end-1)];
        [FEn,M2] = f_Sz2(FE,M1,M2,S1_E,S2_E,N1,N2);
        
        %更新滤波器系数
        W1 = W1 - miu1 * FEn * X2;
       
    end
    E_main = E_main + e_main;
end

E_main = E_main /T;
figure(1);
set(gcf,'Position',[100 100 1200 800])
plot(10*log10(E_main(1:end)),'k-','LineWidth',1.5);
hold on
set(gca,'Ylim',[-SNR-5 5]);



%%  FE_DVM
set(gca,'Ylim',[-40 5]);
for miu3 =2e-3
    
    E_main = zeros(M,1);
    S1_E = flipud(S1);
    S2_E = flipud(S2);
    Beita = N1+N2-2;
    miu1 = 3e-2;
    miu2 = 3E-3;
    yita = 0.02;
    yita2 = 0.02;
    
    for k = 1:T
        
        x_index = 1:M;
        xx = x(x_index);
        d = [dd(x_index(1:M/2));-dd(x_index(M/2+1:M))];
        d = awgn(d,SNR);
        
        x_input = zeros(N+N,1);
        
        X1n = zeros(1,1);
        X1 = zeros(N,1);
        W1 = zeros(N,1);
        
        X2n = zeros(N,1);
        W2 = zeros(N,2);
        W2(1,1) = 1;
        W2(1,2) = 0;
        
        W3 = zeros(N,3);
        W3(1,1) = 1;
        W3(1,2) = 0.5;
        W3(1,3) = 0;
        
        Y2 = zeros(N+1,1);
        Y3 = zeros(N+1,1);
        
        Yn = 0;
        Y = zeros(N,1);
        FFE = zeros(N1+N2,1);
        M1 = zeros(N1+N2,1);
        M2 = zeros(N1+N2,1);
        
        e = zeros(M,1);
        e_main = zeros(M,1);
        ee = EEEd ;
        ddd = EEEd ;
        
        A = zeros(2,1);
        B = zeros(2,1);
        y = zeros(2,1);
        ys = zeros(2,1);
        yos = zeros(2,1);
        ys1 = zeros(N+1,1);
        ys2 = zeros(N+1,1);
        
        A2 = zeros(3,1);
        B2 = zeros(3,1);
        y2 = zeros(3,1);
        yss2 = zeros(3,1);
        yos2 = zeros(3,1);
        ys21 = zeros(N+1,1);
        ys22 = zeros(N+1,1);
        ys23 = zeros(N+1,1);
        u = zeros(2*N,1);
        energy_y = zeros(M,1);
        
        
        i = 1;
        for i = 1:M
            
            %产生x_input 信号
            x_input = [x(i);x_input(1:end-1)];
            
            %产生Volterra X信号
            X1 = x_input(1:N);
            X2 = x_input(Beita+1:Beita+N);
            
            %产生输出信号
            Yn = X1' * W1;
            Yn2 = 1;
            for a = 1:2
                y(a) = X1' * W2(:,a);
                Yn2 = Yn2 * X1' * W2(:,a);
            end
            Yn3 = 1;
            for a = 1:3
                y(a) = X1' * W3(:,a);
                Yn3 = Yn3 * X1' * W3(:,a);
            end
            
            Yn = Yn+Yn2+Yn3;
            Y = [Yn;Y(1:end-1)];
            Y2 = [Yn2;Y(1:end-1)];
            Y3 = [Yn3;Y(1:end-1)];
            
            [Y_s,M1] = f_Sz(Y,M1,S1,S2,N1,N2);
            e(i) = d(i) + Y_s;
            
            ee = 0.95 * ee + 0.05 * e(i)^2;
            ddd = 0.95 * ddd + 0.05 * d(i)^2;
            e_main(i) = ee/ddd;
            
            %产生滤波信号
            FFE = [e(i);FFE(1:end-1)];
            [FEn,M2] = f_Sz2(FFE,M1,M2,S1_E,S2_E,N1,N2);
            
            for a = 1:2
                yos(a) = X1' * W2(:,a);
            end
            A(1) = 1; A(2) = yos(1);
            B(1) = yos(2); B(2) = 1;
            for a = 1:2
                ys(a) = A(a)*B(a);
            end
            
            ys1 = [ys(1);ys1(1:end-1)];
            ys2 = [ys(2);ys2(1:end-1)];
                       
            for a = 1:3
                yos2(a) = X1' * W3(:,a);
            end
            
            A2(1) = 1; A2(2) = yos2(1); A2(3) = yos2(1)*yos2(2);
            B2(1) = yos2(3)*yos2(2); B2(2) = yos2(3);  B2(3) = 1;
            
            for a = 1:3
                yss2(a) = A2(a)*B2(a);
            end
            
            ys21 = [yss2(1);ys21(1:end-1)];
            ys22 = [yss2(2);ys22(1:end-1)];
            ys23 = [yss2(3);ys23(1:end-1)];
            
            %更新滤波器系数
            W1 = W1 - miu1 * FEn * X2;
            
            W2(:,1) = W2(:,1) - miu2 * FEn * X2 * ys1(Beita+1)/(Y2(1)^2*N + yita);
            W2(:,2) = W2(:,2) - miu2 * FEn * X2 * ys2(Beita+1)/(Y2(1)^2*N  + yita);

            W3(:,1) = W3(:,1) - miu3 * FEn * X2 * ys21(Beita+1)/(Y3(1)^2*N  + yita2);
            W3(:,2) = W3(:,2) - miu3 * FEn * X2 * ys22(Beita+1)/(Y3(1)^2*N  + yita2);
            W3(:,3) = W3(:,3) - miu3 * FEn * X2 * ys23(Beita+1)/(Y3(1)^2*N  + yita2);
            
        end
        
        E_main = E_main + e_main;
        
    end
    
    figure(1);
    E_main = E_main /T;
    plot(10*log10(E_main(1:end)),'r-d','LineWidth',1.5,'MarkerIndices',1000:2000:length(E_main),'MarkerSize',8);
    hold on;

    
end


grid on;
xlabel('迭代次数');
ylabel('MSE/dB');
grid on;









