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
xi = data(:, 9); % Lymphocyte Count
yi = data(:, 10); % Neutrophils Count

% Plota um gráfico de dispersão dos dados (pontos)
plot(xi, yi, 'o')
xlim([0, max(xi) * 1.1]) % ajuste do limite no eixo x (margem de 10%)
ylim([0, max(yi) * 1.1]) % ajuste do limite no eixo y (margem de 10%)
grid on
xlabel('Quantidade de Linfoócito')
ylabel('Quantidade de Neutrófilo')

% Cálculo dos parâmetros da regressão linear
% Para encontrar a reta de regressão y=a0 + a1x, precisamos determinar os coeficientes a0 (intercepto) e a1 (inclinação)
n = length(xi); % Número de pontos de dados (número de pacientes)

% Cálculo de a1 aplicando o método de mínimos quadrados, que minimiza a soma dos erros ao quadrado entre os pontos e a reta de ajuste
a1 = (n * sum(xi .* yi) - sum(xi) * sum(yi)) / (n * sum(xi .^ 2) - (sum(xi) ^ 2));
% Cálculo de a0 garantindo que a linha de regressão é centralizada nos dados e minimiza a distância total entre a linha e os pontos
a0 = mean(yi) - a1 * mean(xi);

% Exibir a reta de regressão
hold on % Mantém o gráfico de dispersão existente para sobrepor a reta de regressão
plot(xi, a1 * xi + a0, 'r') # Plota a reta de regressão em vermelho ('r')
title('Regressão Linear de qunatidade de linfoócito versus quantidade de neutroófilo')

% Cálculo dos erros e coeficiente de determinação
St = sum((yi - mean(yi)) .^ 2);  % Soma total dos quadrados
Sr = sum((yi - (a0 + a1 * xi)) .^ 2);  % Soma dos quadrados dos resíduos
r2 = (St - Sr) / St;  % Coeficiente de determinação R²
s_yx = sqrt(Sr / (n - 2));  % Erro padrão da estimativa
s_y = sqrt(St / (n - 1));  % Desvio padrão de yi

% Exibir resultados
fprintf('Coeficientes da regressão: a0 = %.4f, a1 = %.4f\n', a0, a1);
fprintf('Coeficiente de determinação R²: %.4f\n', r2);
fprintf('Erro padrão da estimativa (s_yx): %.4f\n', s_yx);
fprintf('Desvio padrão de yi (s_y): %.4f\n', s_y);


