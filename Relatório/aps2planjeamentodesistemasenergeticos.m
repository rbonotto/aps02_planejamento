%Planejamento de Sistemas Energ�ticos - Prof� G�remi G. Dranka
%APS 2 - FLuxo de Pot�ncia
%Discentes: Rafael Bonotto 
%           Raphael Henrique Soares Machado

close all;
clear all;
clc;
format shortEng
format compact

%% Parte 1 Fluxo de Pot�ncia Linearizado

%Vari�veis do sistema
P1 = 0;
P2 = 0;
P3 = 0;
Pgen1 = 0;
Pgen2 = 0;
Pgen3 = 0;
teta1 = 0;
teta2 = 0;
teta3 = 0;
lambda1 = 0;
lambda2 = 0;
lambda3 = 0;
lambda4 = 0;
lambda5 = 0;

%Reat�ncia da linha por unidade

x12 = 0.1;
x13= 0.125;
x23 = 0.2;

%Declara��o das Equa��es

custo1 = 561 + 7.92*P1 + 0.001562*P1^2;
custo2 = 310 + 7.85*P2 + 0.00194*P2^2;
custo3 = 078 + 7.97*P3 + 0.00482*P3^2;

%Matriz Adimit�ncia

B11= 1/x12 + 1/x13;
B12 = - 1/x12;
B13 = - 1/x13;
B21 = B12;
B22 = 1/x12 + 1/x23;
B23 = - 1/x23;
B31 = B13;
B32 = B23;
B33 = 1/x13 + 1/x23;

B = [B11 B12 B13; B21 B22 B23; B31 B32 B33];

%Resposta

A1 = [0.003124 0 0; 0 .00388 0; 0 0 .00964];
A2 = [1800 -1000 -800; -1000 1500 -500; -800 -500 1300];
A3 = [ 1 0; 0 0; 0 0];
% A3 = [ 1 -B12; 0 B12; 0 0];
A4 =  [0 1; 0 -1; 0 0];

% A = [A1 zeros(3) -eye(3) zeros(3,1); zeros(3) zeros(3) A2 A3; -eye(3) A2 zeros(3) zeros(3,1); zeros(1,3) A3' zeros(1,3) 0]
A = [A1 zeros(3) -eye(3) A4; zeros(3) zeros(3) A2 A3; -eye(3) A2 zeros(3) zeros(3,2); (A4)' (A3)' zeros(2,3) zeros(2)]

% B= [-7.92; -7.85; -7.97; 0; 0; 0; -200; -550; -100; 0];
B= [-7.92; -7.85; -7.97; 0; 0; 0; -200; -550; -100; 0; -200];

resposta = inv(A)*B

Pg1 = resposta(1)
Pc1 = -B(7)
P1 = Pg1 - Pc1
Pg2 = resposta(2)
Pc2 = -B(8)
P2 = Pg2 - Pc2
Pg3 = resposta(3)
Pc3 = -B(9)
P3 = Pg3 - Pc3

P12 = 100*(resposta(4) - resposta(5))/ x12
P13 = 100*(resposta(4) - resposta(6))/ x13
P23 = 100*(resposta(5) - resposta(6))/ x23

%% Parte 2 Fluxo de Pot�ncia Linearizado: Condi��o normal e de emerg�ncia

clear all;
clc;

% N�mero de barras
N =  input('Digite o n�mero de Barras: ');barra_ref = input('Digite o n�mero da barra de refer�ncia: ');

%Declara��o de vari�veis base

%Matriz de reat�ncias
X= zeros(N);
%Matriz Suscept�ncia
B= zeros(N);
%Matriz do fluxo de pot�ncia das linhas
P_linha = zeros(N);
%Vetor das Pot�ncias Geradas nas barras
P_GER = zeros(1,N);
%Vetor das Cargas consumidas nas barras
P_CAR = zeros(1,N);
%Vetor de Pot�ncia l�quida das barras
P = zeros(1,N);
%Vetor de �ngulo das barras
angulo_barra = 1:N;
teta_nome = zeros(N);
% La�o para os valores de entrada das reat�ncias no sistema
for i= 1:N
    for j=(i+1):N
        % if-else para gerar os termos fora da diagonal principal
        if i~=j
            fprintf('Insira a reat�ncia entra a linha %d - %d: ', i, j)
            X(i,j)= input('');
        else
            X(i,i)= 0;
        end
    end
end
% La�o para os valores de entrada das pot�ncias e cargas no sistema
for i= 1:N
    fprintf('Insira o valor da gera��o na barra %d: ', i)
    P_GER(i)= input('');
    fprintf('Insira o valor da carga na barra %d: ', i)
    P_CAR(i)= input('');
end
%Vari�vel de controle para deixar o programa em loop possibilitando a
%escolha de mais cen�rios pelo usu�rio
novo_cenario=1;
while (novo_cenario==1)
    disp('-----------------Escolha o Cen�rio---------------')
    disp('[1] - Considere o caso base (sem conting�ncia)')
    disp('[2] - Conting�ncia em cada linha (N � 1)')
    disp('[3] - Conting�ncia em cada linha (N � 1) com carga extra(%)')
    cenario= input('Escolha o cen�rio: ');
    if cenario==1 %Primeiro Cen�rio
        %Calculando a pot�ncia injetada nas barrras
        P = P_GER - P_CAR;
        %Retirando a barra de refer�ncia dos vetores
        P_inj = P;
        P_inj(barra_ref) = [];
        %Matriz de reat�ncias
        X = X + X';
        % La�o para a forma��o da matriz suscept�ncia no sistema
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
        %Gerando a Matriz B_linha(sem a barra de refer�ncia)
        B_linha=B;
        %Eliminando a linha referente a barra de refer�ncia
        B_linha(barra_ref, :) = [];
        %Eliminando a coluna referente a barra de refer�ncia
        B_linha(:,barra_ref) = [];
        % C�lculo dos �ngulos da barras
        Teta = B_linha\(P_inj)';
        %Vari�vel Teta1 utilizada para inserir o �ngulo de referencia no 
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
        %Visualiza��o dos �ngulos de barra
        Teta_resultados = [ (angulo_barra)' Teta];
        disp('======Teta=======Valores======')
        disp(Teta_resultados)
        %La�o para o c�lculo do fluxo de pot�ncia das linhas
        for ii= 1:N
            for jj= (ii+1):N
                if (X(ii,jj)~=0)
                    fprintf('Fluxo de Pot�ncia na Linha %d-%d(sem conting�ncia): ', ii, jj)
                    P_linha(ii,jj)= (Teta(ii) - Teta(jj))/X(ii,jj);
                    disp(P_linha(ii,jj));
                end
            end
        end
        clear B;
    elseif cenario==2 %Segundo Cen�rio
        %Calculando a pot�ncia injetada nas barrras
        P = P_GER - P_CAR;
        %Retirando a barra de refer�ncia dos vetores
        P_inj = P;
        P_inj(barra_ref) = [];
        %Gerando a matriz de reatancia-Contig�ncia
        X_ctg = X;
        %La�o para a elimina��o de linhas do sistema
        for i= 1:N
            
            for j=((i+1):N)
                % La�o que verifica se a reatancia � diferente de zero
                if ( X_ctg(i,j) ~= 0)
                    fprintf('Linha %d - %d em Falta  \n', i, j);
                    X_ctg(i,j) = inf;
                    X_ctg = X_ctg + (X_ctg)';
                    
                    % La�o para a forma��o da matriz suscept�ncia no
                    % sistema com a retirada da Linha i-j
                    for ii= 1:N
                        for jj=1:N
                            % La�o para termos fora da diagonal principal
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
                    %Gerando a Matriz B_linha(sem a barra de refer�ncia)
                    B_linha=B;
                    %Eliminando a linha referente a barra de refer�ncia
                    B_linha(barra_ref, :) = [];
                    %Eliminando a coluna referente a barra de refer�ncia
                    B_linha(:,barra_ref) = [];
                    % C�lculo dos �ngulos da barras
                    Teta = (B_linha)\(P_inj)';
                    %Vari�vel Teta1 utilizada para inserir o �ngulo de 
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
                    %Visualiza��o dos �ngulos de barra
                    Teta_resultados = [ (angulo_barra)' Teta];
                    disp('======Teta=======Valores======')
                    disp(Teta_resultados)
                    %La�o para o c�lculo do fluxo de pot�ncia das linhas
                    for ii= 1:N
                        for jj= (ii+1):N
                            if (X_ctg(ii,jj)~=0)
                                fprintf('Fluxo da Linha %d-%d com conting�ncia da Linha %d-%d: ', ii, jj,i,j);
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
        
    else %Terceiro Cen�rio
        
        carga_extra= input('Escolha o valor de carga extra(em porcentagem): ');
        carga_extra = 1 + carga_extra/100;
        %Novos valores de carga nas barras
        P_CAR1 = P_CAR * carga_extra;
        %Calculando a pot�ncia injetada nas barrras
        P = P_GER - P_CAR1;
        %Retirando a barra de refer�ncia dos vetores
        P_inj = P;
        P_inj(barra_ref) = [];
        %Gerando a matriz de reatancia-Contig�ncia
        X_ctg = X;
        %La�o para a entrada de reat�ncias no sistema
        for i= 1:N
            for j=((i+1):N)
                if ( X(i,j) ~= 0)
                    fprintf('Linha de Falta %d - %d \n', i, j);
                    X_ctg(i,j) = inf;
                    X_ctg = X_ctg + (X_ctg)';
                    % La�o para a forma��o da matriz suscept�ncia do sistema
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
                    
                    %Gerando a Matriz B_linha(sem a barra de refer�ncia)
                    B_linha=B;
                    %Eliminando a linha referente a barra de refer�ncia
                    B_linha(barra_ref, :) = [];
                    %Eliminando a coluna referente a barra de refer�ncia
                    B_linha(:,barra_ref) = [];
                    % C�lculo dos �ngulos da barras
                    Teta = (B_linha)\(P_inj)';
                    %Vari�vel Teta1 utilizada para inserir o �ngulo de 
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
                    %Visualiza��o dos �ngulos de barra
                    Teta_resultados = [ (angulo_barra)' Teta];
                    disp('======Teta=======Valores======')
                    disp(Teta_resultados)
                    %La�o para o c�lculo do fluxo de pot�ncia das linhas
                    for ii= 1:N
                        for jj= (ii+1):N
                            if (X_ctg(ii,jj)~=0)
                                fprintf('Fluxo de Pot�ncia da Linha %d-%d com conting�ncia da Linha %d-%d: ', ii, jj,i,j);
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
    novo_cenario = input('Deseja escolher um novo cen�rio? [1] - Sim, [0] - N�o: ');
    %Reseta vari�veis base para as equa��es
    if novo_cenario==1
        clear X_ctg;
        clear Teta;
        clear P;
        clear B;
    end
end
