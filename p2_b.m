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
### 1.2

clear all
close all
clc

% Nomes das colunas para referência
columns_names = {'Outcome', 'Patient Age', 'Gender', ...
                 'Ventilated (Y/N)', 'Red blood cell distribution width', ...
                 'Monocytes(%)', 'White blood cell count', ...
                 'Platelet Count', 'Lymphocyte Count', ...
                 'Neutrophils Count', 'Days Hospitalized'};

% Leitura dos dados
data = csvread('COVID-19_CBC_Data_cleaned.csv');

% Removendo a primeira linha que contém strings dos títulos
data(1, :) = [];

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
xlim([0, max(xi_1) * 1.1]) % Ajuste do limite no eixo x
ylim([0, max(yi_1) * 1.1]) % Ajuste do limite no eixo y
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
xlim([0, max(xi_2) * 1.1]) % Ajuste do limite no eixo x
ylim([0, max(yi_2) * 1.1]) % Ajuste do limite no eixo y
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

% Exibir resultados
fprintf('Resultados das análises numéricas pra Regressão Linear de Linfócitos vs Neutrófilos\n')
fprintf('Coeficientes da regressão: a0 = %.4f, a1 = %.4f\n', a0_1, a1_1);
fprintf('Coeficiente de determinação R²: %.4f\n', r2_1);
fprintf('Erro padrão da estimativa (s_yx): %.4f\n', s_yx_1);
fprintf('Desvio padrão de yi (s_y): %.4f\n', s_y_1);
if s_yx_1 < s_y_1
  fprintf('O modelo de regressão é bom!\n')
  else
    fprint('O modelo de regressão não é bom\n')
end

fprintf('Resultados das análises numéricas pra Regressão Linear de Plaqueta vs Dias Hospitalizado\n')
fprintf('Coeficientes da regressão: a0 = %.4f, a1 = %.4f\n', a0_2, a1_2);
fprintf('Coeficiente de determinação R²: %.4f\n', r2_2);
fprintf('Erro padrão da estimativa (s_yx): %.4f\n', s_yx_2);
fprintf('Desvio padrão de yi (s_y): %.4f\n', s_y_2);
if s_yx_2 < s_y_2
  fprintf('O modelo de regressão é bom!\n')
  else
    fprint('O modelo de regressão não é bom\n')
end

fprintf('Comparação entre os dois modelos de regressão\n')
if s_yx_1 > s_yx_2
  fprintf('O modelo de regressão de Plaqueta vs Dias Hospitalizado é melhor que o moedelo de regressão de Linfócitos vs Neutrófilos!\n')
  else
    fprintf('O modelo de regressão de Linfócitos vs Neutrófilos é melhor que o moedelo de regressão de Plaqueta vs Dias Hospitalizado!\n')
 end
