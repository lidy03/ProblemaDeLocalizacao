library(lpSolve)
library(geobr)
library(sf)
library(ggplot2)
library(dplyr)
library(units)

# Baixar e preparar os dados dos bairros de Itacoatiara
bairros <- read_neighborhood(year = 2022)
bairros_itacoatiara <- bairros %>%
  filter(code_muni == 1301902) %>%
  select(code_neighborhood, name_neighborhood, geom) %>%
  mutate(label_id = row_number())

if (nrow(bairros_itacoatiara) == 0) {
  stop("Não foram encontrados bairros para Itacoatiara-AM")
}

# Calcular as coordenadas dos centroides dos bairros e a matriz de distâncias
coords <- st_centroid(bairros_itacoatiara$geom)
coords_matrix <- st_coordinates(coords)

# Verificar se coords_matrix é válida antes de calcular a distância
if (any(is.na(coords_matrix)) || any(is.infinite(coords_matrix))) {
  stop("Coordenadas dos centroides contêm valores inválidos (NA/Inf). Verifique 'st_centroid'.")
}

dist_matrix <- as.matrix(dist(coords_matrix))

# Parâmetros do problema p-centro
p <- 7 # Número de centros a serem selecionados
n <- nrow(bairros_itacoatiara) # Número de bairros (demandas e potenciais locais para centros)

# Verificar se n é válido
if (n == 0) {
  stop("Número de bairros é zero. Não é possível executar o modelo p-centro.")
}

# Implementar o modelo p-centro

# Variáveis de decisão: y_ij (n*n), x_j (n), w (1)
# y_ij: 1 se demanda i é atendida pelo centro j, 0 caso contrário. (n*n variáveis)
# x_j: 1 se um centro é aberto no bairro j, 0 caso contrário. (n variáveis)
# w: a distância máxima de cobertura. (1 variável)

num_vars <- n * n + n + 1

# Função objetivo: minimizar w (a última variável)
obj <- c(rep(0, n * n), # Coeficientes para y_ij
         rep(0, n),     # Coeficientes para x_j
         1)             # Coeficiente para w

# ----------------------------------------------------------------------
# Construção das restrições
# ----------------------------------------------------------------------

# 1. Restrição: Cada demanda i deve ser atendida por exatamente um centro j (Sum_j y_ij = 1)
# n restrições
A1 <- matrix(0, nrow = n, ncol = num_vars)
b1 <- rep(1, n)
dir1 <- rep("=", n)
for (i in 1:n) {
  A1[i, ((i - 1) * n + 1):(i * n)] <- 1 # Y_i1, Y_i2, ..., Y_in
}

# 2. Restrição: Exatamente 'p' centros devem ser selecionados (Sum_j x_j = p)
# 1 restrição
A2 <- matrix(0, nrow = 1, ncol = num_vars)
b2 <- p
dir2 <- "="
A2[1, (n * n + 1):(n * n + n)] <- 1 # X_1, X_2, ..., X_n

# 3. Restrição: y_ij <= x_j (Se i é atendido por j, então j deve ser um centro)
# n*n restrições
A3 <- matrix(0, nrow = n * n, ncol = num_vars)
b3 <- rep(0, n * n)
dir3 <- rep("<=", n * n)
for (i in 1:n) {
  for (j in 1:n) {
    row_idx <- (i - 1) * n + j # Converte (i,j) para um índice linear de 1 a n*n
    col_idx_y_ij <- (i - 1) * n + j # Posição da variável y_ij
    col_idx_x_j <- n * n + j # Posição da variável x_j
    
    A3[row_idx, col_idx_y_ij] <- 1
    A3[row_idx, col_idx_x_j] <- -1
  }
}

# 4. Restrição: d_ij * y_ij <= w para todo i, j (a distância para cada demanda atendida não pode exceder w)
# Nota: Esta restrição deve ser reformulada para o p-centro padrão. A formulação padrão é:
# w >= d_ik * Y_ik para cada k centro de um bairro i.
# Ou seja, a distância da demanda i ao centro que a atende (j) deve ser menor ou igual a w.
# A formulação mais comum para p-centro é "Min w s.t. for each i, w >= min_j(d_ij if X_j=1)"
# Que é equivalente a: Para cada i, w >= d_ij * Y_ij para todo j.
# A formulação que você colocou: w >= Sum_j (d_ij * Y_ij) é para p-mediana, ou outra variante.
# O objetivo do p-centro é garantir que NENHUMA demanda esteja a uma distância maior que 'w'.
# Portanto, para cada demanda 'i', a distância dela ao seu centro alocado (Y_ij=1) deve ser <= w.
# Recomendo a seguinte restrição para p-centro, para cada demanda 'i':
# d_ij * Y_ij - w <= 0 para TODO par (i,j)
# Se Y_ij = 0, a restrição 0 - w <= 0 (ou -w <= 0, w >= 0)
# Se Y_ij = 1, a restrição d_ij - w <= 0 (ou d_ij <= w)
# Isso garante que w é maior ou igual à distância para *qualquer* par i,j que está conectado.
# O 'min' da função objetivo garantirá que w seja a MAIOR dessas distâncias.

A4 <- matrix(0, nrow = n * n, ncol = num_vars) # n*n restrições
b4 <- rep(0, n * n)
dir4 <- rep("<=", n * n)
for (i in 1:n) {
  for (j in 1:n) {
    row_idx <- (i - 1) * n + j
    col_idx_y_ij <- (i - 1) * n + j
    col_idx_w <- num_vars # w é a última variável
    
    A4[row_idx, col_idx_y_ij] <- dist_matrix[i, j] # Coeficiente para d_ij * y_ij
    A4[row_idx, col_idx_w] <- -1                   # Coeficiente para -w
  }
}

# Combinar todas as matrizes de restrição e vetores
A <- rbind(A1, A2, A3, A4)
b <- c(b1, b2, b3, b4)
dir <- c(dir1, dir2, dir3, dir4)

# Verificar se as dimensões estão corretas
if (ncol(A) != num_vars || nrow(A) != length(b) || nrow(A) != length(dir)) {
  stop("Erro de dimensão nas matrizes de restrição. Verifique a construção de A, b, dir.")
}


# Resolver o problema de otimização
# indices_bin são as variáveis y_ij (n*n) e x_j (n)
indices_bin <- 1:(n * n + n)

# A variável w (d_max) não deve ser binária, mas sim contínua.
# lpSolve por padrão trata as variáveis não listadas em int.vec como contínuas.
resultado_p_centro <- lp("min", obj, A, dir, b, int.vec = indices_bin)


# Processar e mapear os resultados do p-centro
bairros_itacoatiara <- bairros_itacoatiara %>%
  mutate(p_centro_selected = "Não Selecionado")

if (resultado_p_centro$status == 0) {
  cat("Solução p-centro encontrada!\n")
  solucao_p_centro <- resultado_p_centro$solution
  
  # Variável renomeada: w_otimo para distancia_max_cobertura_otima
  # w é a última variável
  distancia_max_cobertura_otima <- solucao_p_centro[num_vars]
  
  cat("Distância máxima de cobertura ótima:", distancia_max_cobertura_otima, "\n")
  
  # As variáveis x_j (seleção de centro) estão nas posições n*n + 1 até n*n + n
  postos_otimos_indices_centro <- which(solucao_p_centro[(n * n + 1):(n * n + n)] == 1)
  
  if (length(postos_otimos_indices_centro) > 0) {
    bairros_otimos_centro <- bairros_itacoatiara$name_neighborhood[postos_otimos_indices_centro]
    cat("Os locais ótimos (p-centro) são nos bairros:\n", paste(bairros_otimos_centro, collapse = ", "), "\n")
    bairros_itacoatiara$p_centro_selected[postos_otimos_indices_centro] <- "Selecionado (Ótimo)"
  } else {
    cat("Nenhum posto foi selecionado pelo modelo p-centro (pode indicar problema ou p=0).\n")
  }
  
} else {
  cat("Não foi possível encontrar uma solução ótima para o p-centro. Status:", resultado_p_centro$status, "\n")
  if (resultado_p_centro$status == 2) {
    cat("O problema é infactível (não há solução). Verifique as restrições e o valor de 'p'.\n")
  } else if (resultado_p_centro$status == 3) {
    cat("O problema é ilimitado (solução infinita). Verifique a função objetivo.\n")
  }
}

# ----------------------------------------------------------------------
# Processamento e visualização de UBS existentes
# ----------------------------------------------------------------------
# Carregar e filtrar dados das UBS existentes
servicos_saude_brasil <- read_health_facilities() %>%
  st_as_sf()

ubs_nomes <- c(
  "UNIDADE BASICA DE SAUDE MANOEL MENDES DA SILVA",
  "UNIDADE BASICA DE SAUDE JOSE RESK MAKLOUF",
  "UNIDADE BASICA DE SAUDE NICOLAS EUTHEMES LEKAKIS NETO",
  "UNIDADE BASICA DE SAUDE SANTO ANTONIO",
  "UNIDADE BASICA DE SAUDE PAULO GOMES DA SILVA",
  "UNIDADE BASICA DE SAUDE BERNARDINO DESSIMONI",
  "UNIDADE PRISIONAL DE ITACOATIARA"
)

servicos_saude <- servicos_saude_brasil %>%
  filter(
    abbrev_state == "AM",
    code_muni == 1301902,
    no_fantasia %in% ubs_nomes
  ) %>%
  rename(name = no_fantasia) %>%
  select(name, no_bairro, geom)

if (nrow(servicos_saude) == 0) {
  warning("Nenhuma UBS encontrada com os nomes especificados em Itacoatiara-AM.")
}

# Realizar a junção espacial para associar serviços de saúde aos bairros
servicos_saude_bairros_spatial <- st_join(servicos_saude, bairros_itacoatiara, join = st_within, left = FALSE) %>%
  select(name, label_id)

bairros_com_servicos <- servicos_saude_bairros_spatial %>%
  distinct(label_id) %>%
  pull(label_id)

# Classificar os bairros para mapeamento com as novas legendas
bairros_itacoatiara <- bairros_itacoatiara %>%
  mutate(
    tem_posto_existente = ifelse(label_id %in% bairros_com_servicos, "Existente", "Não Existente"),
    category = case_when(
      p_centro_selected == "Selecionado (Ótimo)" & tem_posto_existente == "Existente" ~ "Atual/Ótima",
      p_centro_selected == "Selecionado (Ótimo)" & tem_posto_existente == "Não Existente" ~ "Localização Ótima",
      tem_posto_existente == "Existente" & p_centro_selected == "Não Selecionado" ~ "Localização UBS Atual",
      TRUE ~ "Nenhum"
    )
  )

# Definir cores para o mapa com as novas legendas
cores_mapa <- c(
  "Atual/Ótima" = "purple",
  "Localização Ótima" = "blue",
  "Localização UBS Atual" = "green",
  "Nenhum" = "lightgrey" # Cor para bairros sem nenhum status
)

# Gerar o mapa
mapa_bairros_numerados_colorido <- ggplot() +
  geom_sf(data = bairros_itacoatiara, aes(fill = category), color = "black", linewidth = 0.2) +
  geom_sf_text(data = bairros_itacoatiara, 
               aes(label = label_id, geometry = st_centroid(geom)),
               size = 2, color = "black", check_overlap = TRUE, fontface = "bold") +
  scale_fill_manual(values = cores_mapa, name = "Legenda") +
  labs(
    title = "Localização de Pontos Ótimos (p-Centro) vs. UBS Existentes em Itacoatiara",
    subtitle = paste0("Distância Máxima Minimizada (w): ", round(distancia_max_cobertura_otima, 2), " unidades")
  ) +
  theme_minimal()

print(mapa_bairros_numerados_colorido)

# Imprimir a legenda numérica dos bairros no console
legenda_numerica_bairros <- bairros_itacoatiara %>%
  st_drop_geometry() %>%
  select(ID = label_id, Bairro = name_neighborhood, Cod_Bairro_IBGE = code_neighborhood) %>%
  arrange(ID)

print("\n--- Legenda dos Bairros ---")
print(legenda_numerica_bairros)
cat("--------------------------\n")