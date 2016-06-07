%Planejamento de Sistemas Energéticos - Profº Géremi G. Dranka
%APS 2 - FLuxo de Potência
%Discentes: Rafael Bonotto 
%           Raphael Henrique Soares Machado

close all;
clear all;
clc;
format shortEng
format compact

%% Parte 1 Fluxo de Potência Linearizado

A = [0.003124 0 0 0 0 0 -1 0 0 0 0; 0 0.00388 0 0 0 0 0 -1 0 0 0; 0 0 2*0.00482 0 0 0 0 0 -1 0 0; 0 0 0 0 0 0 1800 -1000 -800 1 1000;0 0 0 0 0 0 -1000 1500 -500 0 -1000;0 0 0 0 0 0 -800 -500 1300 0 0;-1 0 0 1800 -1000 -800 0 0 0 0 0;0 -1 0 -1000 1500 -500 0 0 0 0 0;0 0 -1 -800 -500 1300 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0 0;0 0 0 1000 -1000 0 0 0 0 0 0];
B = [-7.92;-7.85;-7.97;0;0;0;-200;-550;-100;0;150];

RESP = inv(A)*B;

PG1 = RESP(1)
PG2 = RESP(2)
PG3 = RESP(3)
Theta1 = RESP(4)
Theta2 = RESP(5)
Theta3 = RESP(6)
Lambda1 = RESP(7)
Lambda2 = RESP(8)
Lambda3 = RESP(9)
Lambda4 = RESP(10)
Lambda5 = RESP(11)

%% Parte 2 Fluxo de Potência Linearizado: Condição normal e de emergência

clear all;
clc;

% Número de barras
N =  input('Digite o número de Barras: ');barra_ref = input('Digite o número da barra de referência: ');

%Declaração de variáveis base

%Matriz de reatâncias
X= zeros(N);
%Matriz Susceptância
B= zeros(N);
%Matriz do fluxo de potência das linhas
P_linha = zeros(N);
%Vetor das Potências Geradas nas barras
P_GER = zeros(1,N);
%Vetor das Cargas consumidas nas barras
P_CAR = zeros(1,N);
%Vetor de Potência líquida das barras
P = zeros(1,N);
%Vetor de Ângulo das barras
angulo_barra = 1:N;
teta_nome = zeros(N);
% Laço para os valores de entrada das reatâncias no sistema
for i= 1:N
    for j=(i+1):N
        % if-else para gerar os termos fora da diagonal principal
        if i~=j
            fprintf('Insira a reatância entra a linha %d - %d: ', i, j)
            X(i,j)= input('');
        else
            X(i,i)= 0;
        end
    end
end
% Laço para os valores de entrada das potências e cargas no sistema
for i= 1:N
    fprintf('Insira o valor da geração na barra %d: ', i)
    P_GER(i)= input('');
    fprintf('Insira o valor da carga na barra %d: ', i)
    P_CAR(i)= input('');
end
%Variável de controle para deixar o programa em loop possibilitando a
%escolha de mais cenários pelo usuário
novo_cenario=1;
while (novo_cenario==1)
    disp('-----------------Escolha o Cenário---------------')
    disp('[1] - Considere o caso base (sem contingência)')
    disp('[2] - Contingência em cada linha (N – 1)')
    disp('[3] - Contingência em cada linha (N – 1) com carga extra (%)')
    cenario= input('Escolha o cenário: ');
    if cenario==1 %Primeiro Cenário
        %Calculando a potência injetada nas barrras
        P = P_GER - P_CAR;
        %Retirando a barra de referência dos vetores
        P_inj = P;
        P_inj(barra_ref) = [];
        %Matriz de reatâncias
        X = X + X';
        % Laço para a formação da matriz susceptância no sistema
        for i= 1:N
            for j=1:N
                if i~=j
                    if X(i,j)==0
                        B(i,j)= 0;
                    else
                        B(i,j)= -1/X(i,j);
                        B(i,i)= B(i,i) + 1/X(i,j);
                    end
                end
            end
        end
        %Gerando a Matriz B_linha(sem a barra de referência)
        B_linha=B;
        %Eliminando a linha referente a barra de referência
        B_linha(barra_ref, :) = [];
        %Eliminando a coluna referente a barra de referência
        B_linha(:,barra_ref) = [];
        % Cálculo dos Ângulos da barras
        Teta = B_linha\(P_inj)';
        %Variável Teta1 utilizada para inserir o ângulo de referencia no 
        %no vetor resposta Teta
        Teta1 = Teta';
        switch barra_ref
            case 1
                Teta = [ 0  Teta1];
            case N
                Teta = [Teta1  0];
            otherwise
                Teta = [ Teta1(1:(barra_ref -1)) 0 Teta1((barra_ref +1):N)];
        end
        Teta=(Teta)';
        %Visualização dos ângulos de barra
        Teta_resultados = [ (angulo_barra)' Teta];
        disp('======Teta=======Valores======')
        disp(Teta_resultados)
        %Laço para o cálculo do fluxo de potência das linhas
        for ii= 1:N
            for jj= (ii+1):N
                if (X(ii,jj)~=0)
                    fprintf('Fluxo de Potência na Linha %d-%d(sem contingência): ', ii, jj)
                    P_linha(ii,jj)= (Teta(ii) - Teta(jj))/X(ii,jj);
                    disp(P_linha(ii,jj));
                end
            end
        end
        clear B;
    elseif cenario==2 %Segundo Cenário
        %Calculando a potência injetada nas barrras
        P = P_GER - P_CAR;
        %Retirando a barra de referência dos vetores
        P_inj = P;
        P_inj(barra_ref) = [];
        %Gerando a matriz de reatancia-Contigência
        X_ctg = X;
        %Laço para a eliminação de linhas do sistema
        for i= 1:N
            
            for j=((i+1):N)
                % Laço que verifica se a reatancia é diferente de zero
                if ( X_ctg(i,j) ~= 0)
                    fprintf('Linha %d - %d em Falta  \n', i, j);
                    X_ctg(i,j) = inf;
                    X_ctg = X_ctg + (X_ctg)';
                    
                    % Laço para a formação da matriz susceptância no
                    % sistema com a retirada da Linha i-j
                    for ii= 1:N
                        for jj=1:N
                            % Laço para termos fora da diagonal principal
                            if ii~=jj
                                if X_ctg(ii,jj)==0
                                    B(ii,jj)= 0;
                                else
                                    B(ii,jj)= -1/X_ctg(ii,jj);
                                    B(ii,ii)= B(ii,ii) + 1/X_ctg(ii,jj);
                                end
                                
                            end
                        end
                    end
                    %Gerando a Matriz B_linha(sem a barra de referência)
                    B_linha=B;
                    %Eliminando a linha referente a barra de referência
                    B_linha(barra_ref, :) = [];
                    %Eliminando a coluna referente a barra de referência
                    B_linha(:,barra_ref) = [];
                    % Cálculo dos Ângulos da barras
                    Teta = (B_linha)\(P_inj)';
                    %Variável Teta1 utilizada para inserir o ângulo de 
                    %referencia no vetor resposta Teta
                    Teta1=Teta';
                    switch barra_ref
                        case 1
                            Teta = [ 0  Teta1];
                        case N
                            Teta = [Teta1  0];
                        otherwise
                            Teta = [ Teta1(1:(barra_ref -1)) 0 Teta1((barra_ref +1):N)];
                    end
                    Teta=(Teta)';
                    %Visualização dos ângulos de barra
                    Teta_resultados = [ (angulo_barra)' Teta];
                    disp('======Teta=======Valores======')
                    disp(Teta_resultados)
                    %Laço para o cálculo do fluxo de potência das linhas
                    for ii= 1:N
                        for jj= (ii+1):N
                            if (X_ctg(ii,jj)~=0)
                                fprintf('Fluxo da Linha %d-%d com contingência da Linha %d-%d: ', ii, jj,i,j);
                                P_linha(ii,jj)= (Teta(ii) - Teta(jj))/X_ctg(ii,jj);
                                disp(P_linha(ii,jj));
                            end
                        end
                    end
                    clear B;
                    X_ctg = X;
                end
            end
        end
        
    else %Terceiro Cenário
        
        carga_extra= input('Escolha o valor de carga extra(em porcentagem): ');
        carga_extra = 1 + carga_extra/100;
        %Novos valores de carga nas barras
        P_CAR1 = P_CAR * carga_extra;
        %Calculando a potência injetada nas barrras
        P = P_GER - P_CAR1;
        %Retirando a barra de referência dos vetores
        P_inj = P;
        P_inj(barra_ref) = [];
        %Gerando a matriz de reatancia-Contigência
        X_ctg = X;
        %Laço para a entrada de reatâncias no sistema
        for i= 1:N
            for j=((i+1):N)
                if ( X(i,j) ~= 0)
                    fprintf('Linha de Falta %d - %d \n', i, j);
                    X_ctg(i,j) = inf;
                    X_ctg = X_ctg + (X_ctg)';
                    % Laço para a formação da matriz susceptância do sistema
                    for ii= 1:N
                        for jj=1:N
                            if ii~=jj
                                if X_ctg(ii,jj)==0
                                    B(ii,jj)= 0;
                                else
                                    B(ii,jj)= -1/X_ctg(ii,jj);
                                    B(ii,ii)= B(ii,ii) + 1/X_ctg(ii,jj);
                                end
                            end
                        end
                    end
                    
                    %Gerando a Matriz B_linha(sem a barra de referência)
                    B_linha=B;
                    %Eliminando a linha referente a barra de referência
                    B_linha(barra_ref, :) = [];
                    %Eliminando a coluna referente a barra de referência
                    B_linha(:,barra_ref) = [];
                    % Cálculo dos Ângulos da barras
                    Teta = (B_linha)\(P_inj)';
                    %Variável Teta1 utilizada para inserir o ângulo de 
                    %referencia no vetor resposta Teta
                    Teta1=Teta';
                    switch barra_ref
                        case 1
                            Teta = [ 0  Teta1];
                        case N
                            Teta = [Teta1  0];
                        otherwise
                            Teta = [ Teta1(1:(barra_ref -1)) 0 Teta1((barra_ref +1):N)];
                    end
                    Teta=(Teta)';
                    %Visualização dos ângulos de barra
                    Teta_resultados = [ (angulo_barra)' Teta];
                    disp('======Teta=======Valores======')
                    disp(Teta_resultados)
                    %Laço para o cálculo do fluxo de potência das linhas
                    for ii= 1:N
                        for jj= (ii+1):N
                            if (X_ctg(ii,jj)~=0)
                                fprintf('Fluxo de Potência da Linha %d-%d com contingência da Linha %d-%d: ', ii, jj,i,j);
                                P_linha(ii,jj)= (Teta(ii) - Teta(jj))/X_ctg(ii,jj);
                                disp(P_linha(ii,jj));
                            end
                        end
                    end
                    clear B;
                    X_ctg = X + X';
                end
            end
        end
        
        
    end
    novo_cenario = input('Deseja escolher um novo cenário? [1] - Sim, [0] - Não: ');
    %Reseta variáveis base para as equações
    if novo_cenario==1
        clear X_ctg;
        clear Teta;
        clear P;
        clear B;
    end
end
