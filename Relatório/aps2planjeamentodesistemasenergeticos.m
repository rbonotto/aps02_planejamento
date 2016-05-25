%Planejamento de Sistemas Energ�ticos - Prof� G�remi G. Dranka
%APS 2 - FLuxo de Pot�ncia
%Discentes: Rafael Bonotto 
%           Raphael Henrique Soares Machado

close all;
clear all;
clc;
format shortEng
format compact
%%Parte 1 Fluxo de Pot�ncia Linearizado

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
A3 = [ 1; 0; 0];
% A3 = [ 1 -B12; 0 B12; 0 0];

A = [A1 zeros(3) -eye(3) zeros(3,1); zeros(3) zeros(3) A2 A3; -eye(3) A2 zeros(3) zeros(3,1); zeros(1,3) A3' zeros(1,3) 0]
% A = [A1 zeros(3) -eye(3) zeros(3,2); zeros(3) zeros(3) A2 A3; -eye(3) A2 zeros(3) zeros(3,2); zeros(2,3) A3' zeros(2,3) zeros(2)]

B= [-7.92; -7.85; -7.97; 0; 0; 0; -200; -550; -100; 0];
% B= [-7.92; -7.85; -7.97; 0; 0; 0; -200; -550; -100; 0; -150];

resposta = inv(A)*B

Pg1 = resposta(1);
Pc1 = -B(7);
P1 = Pg1 - Pc1
Pg2 = resposta(2);
Pc2 = -B(8);
P2 = Pg2 - Pc2
Pg3 = resposta(3);
Pc3 = -B(9);
P3 = Pg3 - Pc3

P12 = 100*(resposta(4) - resposta(5))/ x12
P13 = 100*(resposta(4) - resposta(6))/ x13
P23 = 100*(resposta(5) - resposta(6))/ x23