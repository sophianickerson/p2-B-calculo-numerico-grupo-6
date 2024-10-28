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

% Nomes das colunas
columns_names = {'Outcome', 'Patient Age', 'Gender', ...
                 'Ventilated (Y/N)', 'Red blood cell distribution width', ...
                 'Monocytes(%)', 'White blood cell count', ...
                 'Platelet Count', 'Lymphocyte Count', ...
                 'Neutrophils Count', 'Days Hospitalized'};

% Carregar os dados
data = csvread('COVID-19_CBC_Data_cleaned.csv', 1, 0);  % Pulando o cabeçalho (1 linha)

% Separar a coluna do desfecho (Outcome) e as variáveis independentes
outcome = data(:, 1);  % Coluna "Outcome"
variables = data(:, 2:end);  % Todas as variáveis exceto o Outcome
num_vars = size(variables, 2);

% Inicializar vetor de correlações
correlations = zeros(num_vars, 1);

% Calcular a correlação de cada variável com o desfecho
for i = 1:num_vars
    correlations(i) = corr(outcome, variables(:, i));
end

% ---------- REGRESSÃO COM LINFOCITO E NEUTRÓFILOS ----------
figure; % Abre uma nova janela de figura para o primeiro gráfico

xi1 = data(:,9); % Lymphocyte Count
yi1 = data(:,10); % Neutrophils Count

plot(xi1, yi1, 'o')
grid on
xlabel('Linfócito')
ylabel('Neutrófilo')

% Mínimos Quadrados
n = length(xi1);

a11 = (n*sum(xi1.*yi1) - sum(xi1)*sum(yi1))/(n*sum(xi1.^2)-(sum(xi1)^2));
a01 = mean(yi1) - a11*mean(xi1);

hold on
plot(xi1, a11*xi1+a01, 'r')

St1 = sum((yi1 - mean(yi1)).^2);
Sr1 = sum((yi1 - (a01 + a11*xi1)).^2);
r21 = (St1 - Sr1) / St1;
s_yx1 = sqrt(Sr1 / (n - 2));
s_y1 = sqrt(St1 / (n - 1));

fprintf('Soma total dos quadrados teste 1: %d\n', St1);
fprintf('Soma dos quadrados de resíduo teste 1: %d\n', Sr1);
fprintf('Coeficiente de determinação teste 1: %d\n', r21);
fprintf('Erro padrão da estimativa teste 1: %d\n', s_yx1);
fprintf('Desvio padrão de y teste 1: %d\n', s_y1);

% ---------- REGRESSÃO LINEAR PLAQUETAS E DIAS ----------
figure; % Abre uma nova janela de figura para o segundo gráfico

xi2 = data(:,11); % Days Hospitalized
yi2 = data(:,8); % Platelet Count

plot(xi2, yi2, 'o')
grid on
xlabel('Dias Hospitalizados')
ylabel('Plaquetas')

% Mínimos Quadrados
n = length(xi2);

a12 = (n*sum(xi2.*yi2) - sum(xi2)*sum(yi2))/(n*sum(xi2.^2)-(sum(xi2)^2));
a02 = mean(yi2) - a12*mean(xi2);

hold on
plot(xi2, a12*xi2 + a02, 'r')

St2 = sum((yi2 - mean(yi2)).^2);
Sr2 = sum((yi2 - (a02 + a12*xi2)).^2);
r22 = (St2 - Sr2) / St2;
s_yx2 = sqrt(Sr2 / (n - 2));
s_y2 = sqrt(St2 / (n - 1));

fprintf('Soma total dos quadrados teste 2: %d\n', St2);
fprintf('Soma dos quadrados de resíduo teste 2: %d\n', Sr2);
fprintf('Coeficiente de determinação teste 2: %d\n', r22);
fprintf('Erro padrão da estimativa teste 2: %d\n', s_yx2);
fprintf('Desvio padrão de y teste 2: %d\n', s_y2);
