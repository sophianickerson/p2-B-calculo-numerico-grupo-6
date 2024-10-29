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

### 2.1 
# Primeiramente, tinha-se o objetivo de comparar diferentes faixas etárias à quantidade de células sanguíneas do paciente.
# Para isso, comparou-se a média de leucócitos, linfócitos e neutrófilos em 3 faixas etárias (0-40, 40-60 e 60+) anos.
# Para uma melhor visualização, utilizou-se gráficos de barras e lineares.
# segue o código

# REGRESSÃO LINEAR
# Consideramos interessante avaliar a relação entre a quantidade de células sanguíneas por paciente com a taxa de mortalidade
# Para isso foi realizada uma regressão linear, considerando os coeficientes de regressão e o coeficiente de correlação
# Segue o código:

% LEUCÓCITOS X OUTCOME

% Definir as variáveis x e y
x = data(:, 7);       % White blood cell count
y = data(:, 1);       % Outcome

% Cálculo das médias de x e y
x_mean = mean(x);
y_mean = mean(y);

% Cálculo do coeficiente de inclinação (beta1) e do intercepto (beta0) manualmente
numerator_beta1 = sum((x - x_mean) .* (y - y_mean));
denominator_beta1 = sum((x - x_mean).^2);
beta1 = numerator_beta1 / denominator_beta1;

% Intercepto beta0
beta0 = y_mean - beta1 * x_mean;

% Resultados da regressão leucócitos
fprintf("Equação da regressão linear: y = %.4f + %.4f * x\n", beta0, beta1);

% Valores preditos (y_hat) usando a equação de regressão
y_hat = beta0 + beta1 * x;

% Cálculo do R² pelos erros
SSE = sum((y - y_hat).^2);         % Soma dos Erros ao Quadrado
SST = sum((y - y_mean).^2);        % Soma Total dos Quadrados
R2 = 1 - (SSE / SST);              % Coeficiente de determinação (R²)

% Exibir o coeficiente de determinação (R²)
fprintf("Coeficiente de determinação (R²): %.4f\n", R2);


# LINFÓCITOS X OUTCOME

# Definir as variáveis x e y
x = data(:, 9);       # White blood cell count
y = data(:, 1);       # Outcome

# Cálculo das médias de x e y
x_mean = mean(x);
y_mean = mean(y);

# Cálculo do coeficiente de inclinação (beta1) e do intercepto (beta0) manualmente
numerator_beta1 = sum((x - x_mean) .* (y - y_mean));
denominator_beta1 = sum((x - x_mean).^2);
beta1 = numerator_beta1 / denominator_beta1;

# Intercepto beta0
beta0 = y_mean - beta1 * x_mean;

# resultados da regressão linfócitos
fprintf("Equação da regressão linear: y = %.4f + %.4f * x\n", beta0, beta1);

% Valores preditos (y_hat) usando a equação de regressão
y_hat = beta0 + beta1 * x;

% Cálculo do R² pelos erros
SSE = sum((y - y_hat).^2);         % Soma dos Erros ao Quadrado
SST = sum((y - y_mean).^2);        % Soma Total dos Quadrados
R2 = 1 - (SSE / SST);              % Coeficiente de determinação (R²)

% Exibir o coeficiente de determinação (R²)
fprintf("Coeficiente de determinação (R²): %.4f\n", R2);


# NEUTRÓFILOS X OUTCOME

# Definir as variáveis x e y
x = data(:, 10);       # White blood cell count
y = data(:, 1);        # Outcome

# Cálculo das médias de x e y
x_mean = mean(x);
y_mean = mean(y);

# Cálculo do coeficiente de inclinação (beta1) e do intercepto (beta0) manualmente
numerator_beta1 = sum((x - x_mean) .* (y - y_mean));
denominator_beta1 = sum((x - x_mean).^2);
beta1 = numerator_beta1 / denominator_beta1;

# Intercepto beta0
beta0 = y_mean - beta1 * x_mean;

# resultados da regressão neutrófilos
fprintf("Equação da regressão linear: y = %.4f + %.4f * x\n", beta0, beta1);

% Valores preditos (y_hat) usando a equação de regressão
y_hat = beta0 + beta1 * x;

% Cálculo do R² pelos erros
SSE = sum((y - y_hat).^2);         % Soma dos Erros ao Quadrado
SST = sum((y - y_mean).^2);        % Soma Total dos Quadrados
R2 = 1 - (SSE / SST);              % Coeficiente de determinação (R²)

% Exibir o coeficiente de determinação (R²)
fprintf("Coeficiente de determinação (R²): %.4f\n", R2);

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

# CONSIDERAÇÕES FINAIS: Anlisando os coeficientes de correlação da regressão linear e multipla, fica claro que quase não há relação entre a taxa de mortalidade e a quantidade de células sanguíneas por paciente
# OBS: para melhor visualização, gostaria de ter feito gráficos de violino, entretanto encontrei dificuldades de instalar essa função para o octave, logo me limitei a gráficos lineares e de barra.

% Análise 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3.1. Fazendo dois modelos de regressão linear, onde o modelo 𝑦1 = 𝑎0,1 + 𝑎1,1𝑥1 e outro 𝑦2 = 𝑎0,2 + 𝑎1,2𝑥2, onde 𝑥1 e 𝑥2 são variáveis
%que melhor preveem os dias hospitalizados

% Carregando os dados
data = csvread('COVID-19_CBC_Data_cleaned.csv', 1, 0); % Ignorar cabeçalho
dias_hospitalizados = data(:, 11); % Coluna de dias hospitalizados

% Selecionando variáveis independentes (ajustável conforme as variáveis mais correlacionadas)
x1 = data(:, 9); % Contagem de linfócitos
x2 = data(:, 10); % Contagem de neutrófilos

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
    S_y = sqrt(St/(length(y_true) -1))
end

% Calcular métricas para o Modelo 1
[Sr1, r2_1, Sy_x1, S_y1] = calcular_metricas(dias_hospitalizados, y1_pred);

% Calcular métricas para o Modelo 2
[Sr2, r2_2, Sy_x2, S_y2] = calcular_metricas(dias_hospitalizados, y2_pred);

% Exibindo os resultados
fprintf('Modelo 1 (y1 = a0,1 + a1,1 * x1): Sr = %.2f, r^2 = %.2f, Sy/x = %.2f, Sy = %2f\n', Sr1, r2_1, Sy_x1, S_y1);
fprintf('Modelo 2 (y2 = a0,2 + a1,2 * x2): Sr = %.2f, r^2 = %.2f, Sy/x = %.2f, Sy = %2f\n', Sr2, r2_2, Sy_x2, S_y2);

% Comparando as métricas de cada modelo com base em Sy/x e Sy

if Sy_x1 < S_y1
  fprintf('O modelo apresenta boa correlação. (Sy/x < Sy)\n');
else
  fprintf('O modelo não apresenta boa correlação. (Sy/x < Sy)\n');
end


% Comparando os modelos com base no r2
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
x1 = data(:, 9); % Variável independente 1: contagem de linfócitos
x2 = data(:, 10); % Variável independente 2: contagem de neutrófilos

% Calcular os termos necessários para montar o sistema linear
n = length(y);
Sx1 = sum(x1);
Sx2 = sum(x2);
Sy = sum(y);
Sx1x1 = sum(x1 .^ 2);
Sx2x2 = sum(x2 .^ 2);
Sx1x2 = sum(x1 .* x2);
Sx1y = sum(x1 .* y);
Sx2y = sum(x2 .* y);

% Montar a matriz do sistema [X] e o vetor [y]
X = [
    n, Sx1, Sx2;
    Sx1, Sx1x1, Sx1x2;
    Sx2, Sx1x2, Sx2x2
];

Y = [Sy; Sx1y; Sx2y];

% Realizar a eliminação de Gauss para transformar X em uma matriz triangular superior
for i = 1:size(X, 1)
    % Normalizar a linha atual
    X(i, :) = X(i, :) / X(i, i);
    Y(i) = Y(i) / X(i, i);

    % Eliminar as entradas abaixo da diagonal na coluna atual
    for j = i+1:size(X, 1)
        factor = X(j, i);
        X(j, :) = X(j, :) - factor * X(i, :);
        Y(j) = Y(j) - factor * Y(i);
    end
end

% Substituição para trás para resolver o sistema
a = zeros(size(Y));
for i = size(X, 1):-1:1
    a(i) = Y(i) - X(i, i+1:end) * a(i+1:end);
end

% Extrair os coeficientes
a0 = a(1);
a1 = a(2);
a2 = a(3);

% Previsão com o modelo ajustado
y_pred = a0 + a1 * x1 + a2 * x2;

% Calcular Sr, r² e Sy/x para o modelo ajustado
Sr = sum((y - y_pred).^2);
St = sum((y - mean(y)).^2);
r2 = 1 - (Sr / St);
Sy_x = sqrt(Sr / (n - 3));

% Exibir os resultados do modelo ajustado
fprintf('Modelo de regressão múltipla (usando Eliminação de Gauss): Sr = %.2f, r^2 = %.2f, Sy/x = %.2f\n', Sr, r2, Sy_x);

