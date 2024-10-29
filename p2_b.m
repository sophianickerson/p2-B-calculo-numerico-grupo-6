clear all
close all
clc

columns_names = {'Outcome', 'Patient Age', 'Gender', ...
                 'Ventilated (Y/N)', 'Red blood cell distribution width', ...
                 'Monocytes(%)', 'White blood cell count', ...
                 'Platelet Count', 'Lymphocyte Count', ...,
                 'Neutrophils Count', 'Days Hospitalized'};

data = csvread('COVID-19_CBC_Data_cleaned.csv');

# removendo a linha com strings dos títulos
data(1, :) = []

# início dos gráficos de dispersão dois a dois para verificação das relações
# esses gráficos são para vocês terem algum código inicial, não é obrigatório utiliza-los.
num_columns = length(columns_names)-2;

outcome_nao_recuperado = data(data(:, 1) == 0, 1:end);
outcome_recuperado = data(data(:, 1) == 1, 1:end);

for ii = 2:num_columns-1
    figure;
    scatter(outcome_recuperado(:, ii), outcome_recuperado(:, ii+1), 'b', 'x'); % Blue crosses for status 1
    hold on;
    scatter(outcome_nao_recuperado(:, ii), outcome_nao_recuperado(:, ii+1), 'r', 'o'); % Red circles for status 0
    xlabel(columns_names{ii});
    ylabel(columns_names{ii+1});
    title(['Scatter plot of ', columns_names{ii}, ' vs ', columns_names{ii+1}]);
    legend('Não recuperado', 'Recuperado');
    grid on
    hold off;
end

### início do código da prova (fique à vontade para comentar o código de plots dos scatters)
### 1.1
### Essas duas duplas (Linfócitos vs. Neutrófilos e Plaquetas vs. Dias Hospitalizado) foram escolhidas pela observação de uma distribuição com padrões lineares visíveis. Isso é indicativo de que esses pares podem apresentar um ajuste razoável para a regressão linear, com os dados distribuídos de forma a sugerir uma relação direta entre as variáveis. Essa análise inicial ajuda a justificar a escolha antes de partir para o cálculo de regressão linear. Linfócitos vs. Neutrófilos: Na análise gráfica, observamos que há uma relação visualmente mais próxima de linearidade entre as variáveis "Linfócitos" e "Neutrófilos". Esses dois parâmetros sanguíneos estão relacionados entre si na resposta inflamatória e no sistema imunológico, o que é frequentemente esperado em pacientes com doenças infecciosas como a COVID-19. O gráfico de dispersão entre essas variáveis mostra uma tendência de variação conjunta, sugerindo que podem ser bons candidatos para a regressão linear. Plaquetas vs. Dias Hospitalizado: Outro par com uma tendência linear notável é o de "Plaquetas" e "Dias Hospitalizado". A contagem de plaquetas e a duração da internação podem ter uma correlação relacionada à gravidade do quadro clínico. Pacientes com contagens de plaquetas alteradas podem apresentar diferentes tempos de recuperação, refletindo na duração da hospitalização. A relação, apesar de não ser tão forte quanto o primeiro par, ainda sugere uma tendência que pode ser modelada por regressão linear.
### 1.2

% Seleção do par de variáveis "Lymphocyte Count" (coluna 9) e "Neutrophils Count" (coluna 10)
xi_1 = data(:, 9); % Lymphocyte Count
yi_1 = data(:, 10); % Neutrophils Count

% ---------- Gráfico 1: Dispersão dos Dados ----------

figure; % Abre uma nova figura para o gráfico de dispersão
plot(xi_1, yi_1, 'o')
xlim([0, max(xi_1) * 1.1]) % Ajuste do limite no eixo x (margem de 10%)
ylim([0, max(yi_1) * 1.1]) % Ajuste do limite no eixo y (margem de 10%)
grid on
xlabel('Quantidade de Linfócito')
ylabel('Quantidade de Neutrófilo')
title('Dispersão de Linfócitos vs Neutrófilos')

% ---------- Cálculo da Regressão Linear ----------

n_1 = length(xi_1); % Número de pontos de dados (número de pacientes)

% Cálculo de a1 e a0
a1_1 = (n_1 * sum(xi_1 .* yi_1) - sum(xi_1) * sum(yi_1)) / (n_1 * sum(xi_1 .^ 2) - (sum(xi_1) ^ 2));
a0_1 = mean(yi_1) - a1_1 * mean(xi_1);

% ---------- Gráfico 2: Reta de Regressão ----------
plot(xi_1, yi_1, 'o') % Gráfico de dispersão
hold on
plot(xi_1, a1_1 * xi_1 + a0_1, 'r') % Plota a reta de regressão em vermelho
xlim([0, max(xi_1) * 1.1]) % Ajuste do limite no eixo x (margem 10%)
ylim([0, max(yi_1) * 1.1]) % Ajuste do limite no eixo y (margem 10%)
xlabel('Quantidade de Linfócito')
ylabel('Quantidade de Neutrófilo')
title('Regressão Linear de Linfócitos vs Neutrófilos')
grid on
hold off

% ---------- Cálculo dos Erros e Coeficiente de Determinação ----------

St_1 = sum((yi_1 - mean(yi_1)) .^ 2);  % Soma total dos quadrados
Sr_1 = sum((yi_1 - (a0_1 + a1_1 * xi_1)) .^ 2);  % Soma dos quadrados dos resíduos
r2_1 = (St_1 - Sr_1) / St_1;  % Coeficiente de determinação R²
s_yx_1 = sqrt(Sr_1 / (n_1 - 2));  % Erro padrão da estimativa
s_y_1 = sqrt(St_1 / (n_1 - 1));  % Desvio padrão de yi

% Seleção do par de variáveis "Platelet Count" (coluna 8) e "Days Hospitalized" (coluna 11)
xi_2 = data(:, 8); % Platelet Count
yi_2 = data(:, 11); % Days Hospitalized

% ---------- Gráfico 1: Dispersão dos Dados ----------

figure; % Abre uma nova figura para o gráfico de dispersão
plot(xi_2, yi_2, 'o')
xlim([0, max(xi_2) * 1.1]) % Ajuste do limite no eixo x (margem de 10%)
ylim([0, max(yi_2) * 1.1]) % Ajuste do limite no eixo y (margem de 10%)
grid on
xlabel('Quantidade de Plaqueta')
ylabel('Dias Hospitalizado')
title('Dispersão de Plaqueta vs Dias Hospitalizado')

% ---------- Cálculo da Regressão Linear ----------

n_2 = length(xi_2); % Número de pontos de dados (número de pacientes)

% Cálculo de a1 e a0
a1_2 = (n_2 * sum(xi_2 .* yi_2) - sum(xi_2) * sum(yi_2)) / (n_2 * sum(xi_2 .^ 2) - (sum(xi_2) ^ 2));
a0_2 = mean(yi_2) - a1_2 * mean(xi_2);

% ---------- Gráfico 2: Reta de Regressão ----------
plot(xi_2, yi_2, 'o') % Gráfico de dispersão
hold on
plot(xi_2, a1_2 * xi_2 + a0_2, 'r') % Plota a reta de regressão em vermelho
xlim([0, max(xi_2) * 1.1]) % Ajuste do limite no eixo x (margem 10%)
ylim([0, max(yi_2) * 1.1]) % Ajuste do limite no eixo y (margem 10%)
xlabel('Quantidade de Plaqueta')
ylabel('Dias Hospitalizado')
title('Regressão Linear de Plaqueta vs Dias Hospitalizado')
grid on
hold off

% ---------- Cálculo dos Erros e Coeficiente de Determinação ----------

St_2 = sum((yi_2 - mean(yi_2)) .^ 2);  % Soma total dos quadrados
Sr_2 = sum((yi_2 - (a0_2 + a1_2 * xi_2)) .^ 2);  % Soma dos quadrados dos resíduos
r2_2 = (St_2 - Sr_2) / St_2;  % Coeficiente de determinação R²
s_yx_2 = sqrt(Sr_2 / (n_2 - 2));  % Erro padrão da estimativa
s_y_2 = sqrt(St_2 / (n_2 - 1));  % Desvio padrão de yi

% Exibir resultados da análise
fprintf('\n===== ANÁLISE 1: SELEÇÃO DAS VARIÁVEIS =====\n');

fprintf('\n--- 1.1 ---\n');
fprintf('Essas duas duplas (Linfócitos vs. Neutrófilos e Plaquetas vs. Dias Hospitalizado) foram escolhidas pela observação de uma distribuição com padrões lineares visíveis. Isso é indicativo de que esses pares podem apresentar um ajuste razoável para a regressão linear, com os dados distribuídos de forma a sugerir uma relação direta entre as variáveis. Essa análise inicial ajuda a justificar a escolha antes de partir para o cálculo de regressão linear. Linfócitos vs. Neutrófilos: Na análise gráfica, observamos que há uma relação visualmente mais próxima de linearidade entre as variáveis "Linfócitos" e "Neutrófilos". Esses dois parâmetros sanguíneos estão relacionados entre si na resposta inflamatória e no sistema imunológico, o que é frequentemente esperado em pacientes com doenças infecciosas como a COVID-19. O gráfico de dispersão entre essas variáveis mostra uma tendência de variação conjunta, sugerindo que podem ser bons candidatos para a regressão linear. Plaquetas vs. Dias Hospitalizado: Outro par com uma tendência linear notável é o de "Plaquetas" e "Dias Hospitalizado". A contagem de plaquetas e a duração da internação podem ter uma correlação relacionada à gravidade do quadro clínico. Pacientes com contagens de plaquetas alteradas podem apresentar diferentes tempos de recuperação, refletindo na duração da hospitalização. A relação, apesar de não ser tão forte quanto o primeiro par, ainda sugere uma tendência que pode ser modelada por regressão linear.\n')

% Resultados para Linfócitos vs Neutrófilos
fprintf('\n--- Resultados da Regressão Linear: Linfócitos vs Neutrófilos ---\n');
fprintf('  Coeficiente a0               = %10.4f\n', a0_1);
fprintf('  Coeficiente a1               = %10.4f\n', a1_1);
fprintf('  Coeficiente de determinação (R²) = %10.4f\n', r2_1);
fprintf('  Erro padrão da estimativa (s_yx) = %10.4f\n', s_yx_1);
fprintf('  Desvio padrão de yi (s_y)    = %10.4f\n', s_y_1);
if s_yx_1 < s_y_1
    fprintf('  >> O modelo de regressão é bom!\n');
else
    fprintf('  >> O modelo de regressão não é bom\n');
end

% Resultados para Plaquetas vs Dias Hospitalizado
fprintf('\n--- Resultados da Regressão Linear: Plaquetas vs Dias Hospitalizado ---\n');
fprintf('  Coeficiente a0               = %10.4f\n', a0_2);
fprintf('  Coeficiente a1               = %10.4f\n', a1_2);
fprintf('  Coeficiente de determinação (R²) = %10.4f\n', r2_2);
fprintf('  Erro padrão da estimativa (s_yx) = %10.4f\n', s_yx_2);
fprintf('  Desvio padrão de yi (s_y)    = %10.4f\n', s_y_2);
if s_yx_2 < s_y_2
    fprintf('  >> O modelo de regressão é bom!\n');
else
    fprintf('  >> O modelo de regressão não é bom\n');
end

% Comparação entre os modelos
fprintf('\n--- Comparação entre os Modelos ---\n');
if r2_2 > r2_1
  fprintf('>> O modelo de regressão de Plaquetas vs Dias Hospitalizado é melhor que o modelo de Linfócitos vs Neutrófilos!\n')
  else
    fprintf('>> O modelo de regressão de Linfócitos vs Neutrófilos é melhor que o modelo de Plaquetas vs Dias Hospitalizado!\n')
 end
 fprintf('=============================================\n');

 fprintf('\n===== ANÁLISE 2: COMPARAÇÃO ESTATÍSTICA DAS MÉTRICAS =====\n');

### 2.1 
# Primeiramente, tinha-se o objetivo de comparar diferentes faixas etárias à quantidade de células sanguíneas do paciente.
# Para isso, comparou-se a média de leucócitos, linfócitos e neutrófilos em 3 faixas etárias (0-40, 40-60 e 60+) anos.
# Para uma melhor visualização, utilizou-se gráficos de barras e lineares.
# segue o código

# REGRESSÃO POLINOMIAL
# Consideramos interessante avaliar a relação entre a quantidade de células sanguíneas por paciente com a taxa de mortalidade
# Para isso foi realizada uma regressão polinomial de terceiro grau
# Segue o código:

% LEUCÓCITOS VS OUTCOME

xi_3 = data(:, 7); % white blood cell count
yi_3 = data(:, 1); % outcome

% ---------- Gráfico 1: Dispersão dos Dados ----------
figure; % Abre uma nova figura para o gráfico de dispersão
plot(xi_3, yi_3, 'o')
xlim([0, max(xi_3) * 1.1]) % Ajuste do limite no eixo x (margem de 10%)
ylim([0, max(yi_3) * 1.1]) % Ajuste do limite no eixo y (margem de 10%)
grid on
xlabel('Quantidade de Leucócitos')
ylabel('Outcome')
title('Dispersão de Leucócitos vs Outcome')

% ---------- Regressão Polinomial (Grau 3) ----------
n_3 = length(xi_3); % Número de pontos de dados (número de pacientes)

% Montar a matriz de design para polinômio de grau 3: [1, x, x^2, x^3]
X_3 = [ones(n_3, 1), xi_3, xi_3 .^ 2, xi_3 .^ 3];

% Resolver o sistema X * beta = y usando decomposição LU
XtX_3 = X_3' * X_3;
Xty_3 = X_3' * yi_3;
[L, U] = lu(XtX_3);

% Resolve o sistema Lz = Xty
z_3 = L \ Xty_3;

% Resolve o sistema Ubeta = z para obter os coeficientes beta
beta_3 = U \ z_3;

% Coeficientes da regressão polinomial
a0_3 = beta_3(1);
a1_3 = beta_3(2);
a2_3 = beta_3(3);
a3_3 = beta_3(4);

% ---------- Gráfico 2: Curva de Regressão Polinomial ----------
xi_3_sorted = sort(xi_3); % Organizar os valores de x para um gráfico mais suave
yi_3_pred = a0_3 + a1_3 * xi_3_sorted + a2_3 * (xi_3_sorted .^ 2) + a3_3 * (xi_3_sorted .^ 3);

figure; % Abre uma nova figura para o gráfico da regressão
plot(xi_3, yi_3, 'o') % Gráfico de dispersão dos dados originais
hold on
plot(xi_3_sorted, yi_3_pred, 'r') % Plota a curva de regressão em vermelho
xlim([0, max(xi_3) * 1.1]) % Ajuste do limite no eixo x (margem de 10%)
ylim([0, max(yi_3) * 1.1]) % Ajuste do limite no eixo y (margem de 10%)
xlabel('Quantidade de Leucócitos')
ylabel('Outcome')
title('Regressão Polinomial (Grau 3) de Leucócitos vs Outcome')
grid on
hold off

% ---------- Cálculo dos Erros e Coeficiente de Determinação ----------
yi_3_fitted = a0_3 + a1_3 * xi_3 + a2_3 * (xi_3 .^ 2) + a3_3 * (xi_3 .^ 3);
St_3 = sum((yi_3 - mean(yi_3)) .^ 2);  % Soma total dos quadrados
Sr_3 = sum((yi_3 - yi_3_fitted) .^ 2);  % Soma dos quadrados dos resíduos
r2_3 = 1 - (Sr_3 / St_3);  % Coeficiente de determinação R²
s_yx_3 = sqrt(Sr_3 / (n_3 - 4));  % Erro padrão da estimativa (ajustado para grau 3)
s_y_3 = sqrt(St_3 / (n_3 - 1));  % Desvio padrão de yi

% Resultados para Leucócitos vs Outcome (Regressão Polinomial)
fprintf('\n--- Resultados da Regressão Polinomial (Grau 3): Leucócitos vs Outcome ---\n');
fprintf('  Coeficiente a0               = %10.4f\n', a0_3);
fprintf('  Coeficiente a1               = %10.4f\n', a1_3);
fprintf('  Coeficiente a2               = %10.4f\n', a2_3);
fprintf('  Coeficiente a3               = %10.4f\n', a3_3);
fprintf('  Coeficiente de determinação (R²) = %10.4f\n', r2_3);
fprintf('  Erro padrão da estimativa (s_yx) = %10.4f\n', s_yx_3);
fprintf('  Desvio padrão de yi (s_y)    = %10.4f\n', s_y_3);
if s_yx_3 < s_y_3
    fprintf('  >> O modelo de regressão polinomial é bom!\n');
else
    fprintf('  >> O modelo de regressão polinomial não é bom\n');
end

% LINFÓCITOS VS OUTCOME

xi_4 = data(:, 9); % lymphocytes count
yi_4 = data(:, 1); % outcome

% ---------- Gráfico 1: Dispersão dos Dados ----------
figure; % Abre uma nova figura para o gráfico de dispersão
plot(xi_4, yi_4, 'o')
xlim([0, max(xi_4) * 1.1]) % Ajuste do limite no eixo x (margem de 10%)
ylim([0, max(yi_4) * 1.1]) % Ajuste do limite no eixo y (margem de 10%)
grid on
xlabel('Quantidade de Linfócitos')
ylabel('Outcome')
title('Dispersão de Linfócitos vs Outcome')

% ---------- Regressão Polinomial (Grau 3) ----------
n_4 = length(xi_4); % Número de pontos de dados (número de pacientes)

% Montar a matriz de design para polinômio de grau 3: [1, x, x^2, x^3]
X_4 = [ones(n_4, 1), xi_4, xi_4 .^ 2, xi_4 .^ 3];

% Resolver o sistema X * beta = y usando decomposição LU
XtX_4 = X_4' * X_4;
Xty_4 = X_4' * yi_4;
[L, U] = lu(XtX_4);

% Resolve o sistema Lz = Xty
z_4 = L \ Xty_4;

% Resolve o sistema Ubeta = z para obter os coeficientes beta
beta_4 = U \ z_4;

% Coeficientes da regressão polinomial
a0_4 = beta_4(1);
a1_4 = beta_4(2);
a2_4 = beta_4(3);
a3_4 = beta_4(4);

% ---------- Gráfico 2: Curva de Regressão Polinomial ----------
xi_4_sorted = sort(xi_4); % Organizar os valores de x para um gráfico mais suave
yi_4_pred = a0_4 + a1_4 * xi_4_sorted + a2_4 * (xi_4_sorted .^ 2) + a3_4 * (xi_4_sorted .^ 3);

figure; % Abre uma nova figura para o gráfico da regressão
plot(xi_4, yi_4, 'o') % Gráfico de dispersão dos dados originais
hold on
plot(xi_4_sorted, yi_4_pred, 'r') % Plota a curva de regressão em vermelho
xlim([0, max(xi_4) * 1.1]) % Ajuste do limite no eixo x (margem de 10%)
ylim([0, max(yi_4) * 1.1]) % Ajuste do limite no eixo y (margem de 10%)
xlabel('Quantidade de Linfócitos')
ylabel('Outcome')
title('Regressão Polinomial (Grau 3) de Linfócitos vs Outcome')
grid on
hold off

% ---------- Cálculo dos Erros e Coeficiente de Determinação ----------
yi_4_fitted = a0_4 + a1_4 * xi_4 + a2_4 * (xi_4 .^ 2) + a3_4 * (xi_4 .^ 3);
St_4 = sum((yi_4 - mean(yi_4)) .^ 2);  % Soma total dos quadrados
Sr_4 = sum((yi_4 - yi_4_fitted) .^ 2);  % Soma dos quadrados dos resíduos
r2_4 = 1 - (Sr_4 / St_4);  % Coeficiente de determinação R²
s_yx_4 = sqrt(Sr_4 / (n_4 - 4));  % Erro padrão da estimativa (ajustado para grau 3)
s_y_4 = sqrt(St_4 / (n_4 - 1));  % Desvio padrão de yi

% Resultados para Linfócitos vs Outcome (Regressão Polinomial)
fprintf('\n--- Resultados da Regressão Polinomial (Grau 3): Linfócitos vs Outcome ---\n');
fprintf('  Coeficiente a0               = %10.4f\n', a0_4);
fprintf('  Coeficiente a1               = %10.4f\n', a1_4);
fprintf('  Coeficiente a2               = %10.4f\n', a2_4);
fprintf('  Coeficiente a3               = %10.4f\n', a3_4);
fprintf('  Coeficiente de determinação (R²) = %10.4f\n', r2_4);
fprintf('  Erro padrão da estimativa (s_yx) = %10.4f\n', s_yx_4);
fprintf('  Desvio padrão de yi (s_y)    = %10.4f\n', s_y_4);
if s_yx_4 < s_y_4
    fprintf('  >> O modelo de regressão polinomial é bom!\n');
else
    fprintf('  >> O modelo de regressão polinomial não é bom\n');
end


% NEUTRÓFILOS VS OUTCOME

xi_5 = data(:, 10); % neutrophils count
yi_5 = data(:, 1); % outcome

% ---------- Gráfico 1: Dispersão dos Dados ----------
figure; % Abre uma nova figura para o gráfico de dispersão
plot(xi_5, yi_5, 'o')
xlim([0, max(xi_5) * 1.1]) % Ajuste do limite no eixo x (margem de 10%)
ylim([0, max(yi_5) * 1.1]) % Ajuste do limite no eixo y (margem de 10%)
grid on
xlabel('Quantidade de Neutrófilos')
ylabel('Outcome')
title('Dispersão de Neutrófilos vs Outcome')

% ---------- Regressão Polinomial (Grau 3) ----------
n_5 = length(xi_5); % Número de pontos de dados (número de pacientes)

% Montar a matriz de design para polinômio de grau 3: [1, x, x^2, x^3]
X_5 = [ones(n_5, 1), xi_5, xi_5 .^ 2, xi_5 .^ 3];

% Resolver o sistema X * beta = y usando decomposição LU
XtX_5 = X_5' * X_5;
Xty_5 = X_5' * yi_5;
[L, U] = lu(XtX_5);

% Resolve o sistema Lz = Xty
z_5 = L \ Xty_5;

% Resolve o sistema Ubeta = z para obter os coeficientes beta
beta_5 = U \ z_5;

% Coeficientes da regressão polinomial
a0_5 = beta_5(1);
a1_5 = beta_5(2);
a2_5 = beta_5(3);
a3_5 = beta_5(4);

% ---------- Gráfico 2: Curva de Regressão Polinomial ----------
xi_5_sorted = sort(xi_5); % Organizar os valores de x para um gráfico mais suave
yi_5_pred = a0_5 + a1_5 * xi_5_sorted + a2_5 * (xi_5_sorted .^ 2) + a3_5 * (xi_5_sorted .^ 3);

figure; % Abre uma nova figura para o gráfico da regressão
plot(xi_5, yi_5, 'o') % Gráfico de dispersão dos dados originais
hold on
plot(xi_5_sorted, yi_5_pred, 'r') % Plota a curva de regressão em vermelho
xlim([0, max(xi_5) * 1.1]) % Ajuste do limite no eixo x (margem de 10%)
ylim([0, max(yi_5) * 1.1]) % Ajuste do limite no eixo y (margem de 10%)
xlabel('Quantidade de Neutrófilos')
ylabel('Outcome')
title('Regressão Polinomial (Grau 3) de Neutrófilos vs Outcome')
grid on
hold off

% ---------- Cálculo dos Erros e Coeficiente de Determinação ----------
yi_5_fitted = a0_5 + a1_5 * xi_5 + a2_5 * (xi_5 .^ 2) + a3_5 * (xi_5 .^ 3);
St_5 = sum((yi_5 - mean(yi_5)) .^ 2);  % Soma total dos quadrados
Sr_5 = sum((yi_5 - yi_5_fitted) .^ 2);  % Soma dos quadrados dos resíduos
r2_5 = 1 - (Sr_5 / St_5);  % Coeficiente de determinação R²
s_yx_5 = sqrt(Sr_5 / (n_5 - 4));  % Erro padrão da estimativa (ajustado para grau 3)
s_y_5 = sqrt(St_5 / (n_5 - 1));  % Desvio padrão de yi

% Resultados para Neutrófilos vs Outcome (Regressão Polinomial)
fprintf('\n--- Resultados da Regressão Polinomial (Grau 3): Neutrófilos vs Outcome ---\n');
fprintf('  Coeficiente a0               = %10.4f\n', a0_5);
fprintf('  Coeficiente a1               = %10.4f\n', a1_5);
fprintf('  Coeficiente a2               = %10.4f\n', a2_5);
fprintf('  Coeficiente a3               = %10.4f\n', a3_5);
fprintf('  Coeficiente de determinação (R²) = %10.4f\n', r2_5);
fprintf('  Erro padrão da estimativa (s_yx) = %10.4f\n', s_yx_5);
fprintf('  Desvio padrão de yi (s_y)    = %10.4f\n', s_y_5);
if s_yx_5 < s_y_5
    fprintf('  >> O modelo de regressão polinomial é bom!\n');
else
    fprintf('  >> O modelo de regressão polinomial não é bom\n');
end

#2.2

columns_names = {'Outcome', 'Patient Age', 'Gender', ...
                 'Ventilated (Y/N)', 'Red blood cell distribution width', ...
                 'Monocytes(%)', 'White blood cell count', ...
                 'Platelet Count', 'Lymphocyte Count', ...,
                 'Neutrophils Count', 'Days Hospitalized'};

data = csvread('COVID-19_CBC_Data_cleaned.csv');

# removendo a linha com strings dos títulos
data(1, :) = [];

% Definir faixas etárias
age_groups = {[0, 40], [40, 60], [60, Inf]};
group_labels = {'0-40', '40-60', '60+'};

% Índices das colunas
age_column = 2;
wbc_column = 7;
neutrophils_column = 10;
lymphocyte_column = 9;

% HISTOGRAMA DE GLÓBULOS BRANCOS x FAIXA ETÁRIA
figure;
for i = 1:length(age_groups)
    age_range = age_groups{i};
    age_filter = data(:, age_column) >= age_range(1) & data(:, age_column) < age_range(2);
    wbc_data = data(age_filter, wbc_column);
    subplot(3, 1, i);
    hist(wbc_data, 10);  % Ajuste o número de bins conforme necessário
    title(['Contagem de Glóbulos Brancos - Faixa Etária' group_labels{i}]);
    xlabel('Contagem de Glóbulos Brancos');
    ylabel('Frequência');
end

% HISTOGRAMA DE NEUTRÓFILOS x FAIXA ETÁRIA
figure;
for i = 1:length(age_groups)
    age_range = age_groups{i};
    age_filter = data(:, age_column) >= age_range(1) & data(:, age_column) < age_range(2);
    neutrophils_data = data(age_filter, neutrophils_column);
    subplot(3, 1, i);
    hist(neutrophils_data, 10);  % Ajuste o número de bins conforme necessário
    title(['Contagem de Neutrófilos - Faixa Etária' group_labels{i}]);
    xlabel('Contagem de Neutrófilos');
    ylabel('Frequência');
end

% HISTOGRAMA DE LINFÓCITOS x FAIXA ETÁRIA
figure;
for i = 1:length(age_groups)
    age_range = age_groups{i};
    age_filter = data(:, age_column) >= age_range(1) & data(:, age_column) < age_range(2);
    lymphocyte_data = data(age_filter, lymphocyte_column);
    subplot(3, 1, i);
    hist(lymphocyte_data, 10);  % Ajuste o número de bins conforme necessário
    title(['Contagem de Linfócitos - Faixa Etária' group_labels{i}]);
    xlabel('Contagem de Linfócitos');
    ylabel('Frequência');
end


% HISTOGRAMA DE OUTCOME POR FAIXA ETÁRIA

% Índices das colunas relevantes
outcome_column = 1;
age_column = 2;

figure;
for i = 1:length(age_groups)
    age_range = age_groups{i};
    age_filter = data(:, age_column) >= age_range(1) & data(:, age_column) < age_range(2);
    outcome_data = data(age_filter, outcome_column);

    % Contar recuperados (1) e não recuperados (0) na faixa etária
    num_recovered = sum(outcome_data == 1);
    num_not_recovered = sum(outcome_data == 0);

    % Dados para o histograma de barras
    bar_data = [num_not_recovered, num_recovered];
    subplot(3, 1, i);
    bar([0, 1], bar_data);

    title(['Outcome - Faixa Etária ' group_labels{i}]);
    xlabel('Outcome (0: Não Recuperado, 1: Recuperado)');
    ylabel('Frequência');
    xticks([0 1]);
end

# CONSIDERAÇÕES FINAIS: Anlisando os coeficientes de correlação da regressão polinomial e multipla, fica claro que quase não há relação entre a taxa de mortalidade e a quantidade de células sanguíneas por paciente
# OBS: para melhor visualização, gostaria de ter feito gráficos de violino, entretanto encontrei dificuldades de instalar essa função para o octave, logo me limitei a gráficos lineares e de barra.

% Análise 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3.1. Fazendo dois modelos de regressão linear, onde o modelo 𝑦1 = 𝑎0,1 + 𝑎1,1𝑥1 e outro 𝑦2 = 𝑎0,2 + 𝑎1,2𝑥2, onde 𝑥1 e 𝑥2 são variáveis
%que melhor preveem os dias hospitalizados

% Carregando os dados
data = csvread('COVID-19_CBC_Data_cleaned.csv', 1, 0); % Ignorar cabeçalho
dias_hospitalizados = data(:, 11); % Coluna de dias hospitalizados

% Selecionando variáveis independentes (ajustável conforme as variáveis mais correlacionadas)
x1 = data(:, 8); % Contagem de plaqueta
x2 = data(:, 2); % Contagem de neutrófilos

% Função para calcular os coeficientes de um modelo linear simples y = a0 + a1*x
function [a0, a1] = regressao_linear(x, y)
    n = length(x);
    x_mean = mean(x);
    y_mean = mean(y);
    a1 = sum((x - x_mean) .* (y - y_mean)) / sum((x - x_mean).^2);
    a0 = y_mean - a1 * x_mean;
end

% Modelo 1: y1 = a0,1 + a1,1 * x1
[a0_1, a1_1] = regressao_linear(x1, dias_hospitalizados);
y1_pred = a0_1 + a1_1 * x1;

% Modelo 2: y2 = a0,2 + a1,2 * x2
[a0_2, a1_2] = regressao_linear(x2, dias_hospitalizados);
y2_pred = a0_2 + a1_2 * x2;

% Função para calcular Sr, r^2, Sy/x e Sy
function [Sr, r2, Sy_x, S_y] = calcular_metricas(y_true, y_pred)
    Sr = sum((y_true - y_pred).^2);
    St = sum((y_true - mean(y_true)).^2);
    r2 = 1 - (Sr / St);
    Sy_x = sqrt(Sr / (length(y_true) - 2));
    S_y = sqrt(St/(length(y_true) -1));
end

% Calcular métricas para o Modelo 1
[Sr1, r2_1, Sy_x1, S_y1] = calcular_metricas(dias_hospitalizados, y1_pred);

% Calcular métricas para o Modelo 2
[Sr2, r2_2, Sy_x2, S_y2] = calcular_metricas(dias_hospitalizados, y2_pred);

% Exibindo os resultados
fprintf('\n===== ANÁLISE 3: PREDIÇÃO =====\n');
fprintf('\n--- Modelo 1: Plaqueta por Dias Hospitalizado ---\n');
fprintf('(y1 = a0,1 + a1,1 * x1): Sr = %.2f, r^2 = %.2f, Sy/x = %.2f, Sy = %2f\n', Sr1, r2_1, Sy_x1, S_y1);

% Comparando as métricas do modelo 1 com base em Sy/x e Sy

if Sy_x1 < S_y1
  fprintf('O modelo 1 apresenta boa correlação. (Sy/x < Sy)\n');
else
  fprintf('O modelo 1 não apresenta boa correlação. (Sy/x < Sy)\n');
end

fprintf('\n--- Modelo 2: Idade do Paciente por Dias Hospitalizado ---\n');
fprintf('(y2 = a0,2 + a1,2 * x2): Sr = %.2f, r^2 = %.2f, Sy/x = %.2f, Sy = %2f\n', Sr2, r2_2, Sy_x2, S_y2);

% Comparando as métricas do modelo 2 com base em Sy/x e Sy

if Sy_x2 < S_y2
  fprintf('O modelo 2 apresenta boa correlação. (Sy/x < Sy)\n');
else
  fprintf('O modelo 2 não apresenta boa correlação. (Sy/x < Sy)\n');
end

% Comparando os modelos com base no r2
fprintf('\n--- Comparação entre o modelos ---\n');
if r2_1 > r2_2
  fprintf('O Modelo 1 é melhor com base em r^2.\n');
  elseif r2_2 > r2_1
    fprintf('O Modelo 2 é melhor com base em r^2.\n');
  elseif Sy_x1> Sy_x2
    fprintf('O Modelo 2 é melhor com base em Sy/x.\n')
  elseif Sy_x1< Sy_x2
    fprintf('O Modelo 1 é melhor com base em Sy/x.\n')
  else
    fprintf('Ambos os modelos têm desempenho semelhante.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3.2. Implementando um terceiro modelo de regressão linear que combine as duas anteriores

% Carregar os dados
data = csvread('COVID-19_CBC_Data_cleaned.csv', 1, 0); % Ignorar cabeçalho
y = data(:, 11); % Variável dependente: dias hospitalizados
x1 = data(:, 8); % Variável independente 1: contagem de plaquetas
x2 = data(:, 2); % Variável independente 2: idade do paciente

% Regressão linear múltipla usando decomposição LU
% para prever dias hospitalizados com as variáveis independentes x1 e x2.

% Função de regressão linear múltipla usando decomposição LU
function [a0, a1, a2] = regressao_linear_multipla(x1, x2, y)
    % Montando a matriz de design com uma coluna de 1's para o intercepto
    X = [ones(length(x1), 1), x1, x2];

    % Calculando XtX e XtY para o sistema linear
    XtX = X' * X;
    XtY = X' * y;

    % Decomposição LU para resolver XtX * a = XtY
    [L, U] = lu(XtX);

    % Resolvendo Lz = XtY para z usando uma função customizada
    z = resolver_sistema(L, XtY);

    % Resolvendo Ua = z para os coeficientes a
    a = resolver_sistema(U, z);

    % Atribuindo coeficientes a0, a1, e a2
    a0 = a(1);
    a1 = a(2);
    a2 = a(3);
end

% Função customizada para resolver sistemas lineares usando L e U
function x = resolver_sistema(M, b)
    x = zeros(size(b));
    n = length(b);

    % Resolver triangular inferior (L)
    if isdiag(tril(M))
        for i = 1:n
            x(i) = (b(i) - M(i, 1:i-1) * x(1:i-1)) / M(i, i);
        end
    % Resolver triangular superior (U)
    elseif isdiag(triu(M))
        for i = n:-1:1
            x(i) = (b(i) - M(i, i+1:end) * x(i+1:end)) / M(i, i);
        end
    end
end

% Calcular os coeficientes para o modelo de regressão múltipla
[a0_3, a1_3, a2_3] = regressao_linear_multipla(x1, x2, dias_hospitalizados);
y3_pred = a0_3 + a1_3 * x1 + a2_3 * x2;

% Calcular métricas para o modelo de regressão múltipla (Modelo 3)
[Sr3, r2_3, Sy_x3] = calcular_metricas(dias_hospitalizados, y3_pred);

% Exibir resultados do Modelo 3
fprintf('\n--- Modelo 2: Plaqueta + Idade do Paciente por Dias Hospitalizado ---\n');
fprintf('(y3 = a0,3 + a1,3 * x1 + a2,3 * x2): Sr = %.2f, r^2 = %.2f, Sy/x = %.2f\n', Sr3, r2_3, Sy_x3);

% Comparação final dos modelos com base no r²
if r2_3 > max([r2_1, r2_2])
    fprintf('O Modelo 3 é o melhor com base em r^2.\n');
elseif r2_1 > r2_2
    fprintf('O Modelo 1 ainda é melhor com base em r^2.\n');
else
    fprintf('O Modelo 2 ainda é melhor com base em r^2.\n');
end

% 3.3
fprintf('\n--- 3.3 ---\n');
fprintf('\n--- De acordo com o ponto 3.1, o modelo, apesar de ser bom para o exercício, não é suficientemente bom para um modelo preditivo. Isso porque o ajuste não é muito bom, com R2= 0,25. Dessa maneira, não é possível prever a quantidade de dias que o paciente ficará hospitalizado ao olhar a quantidade de plaquetas. ---\n');
