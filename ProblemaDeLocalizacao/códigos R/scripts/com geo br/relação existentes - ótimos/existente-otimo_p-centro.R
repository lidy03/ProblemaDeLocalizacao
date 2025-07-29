library(lpSolve)
library(geobr)
library(sf)
library(ggplot2)
library(dplyr)
library(units)

#Baixar e preparar os dados dos bairros de Itacoatiara
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
dist_matrix <- as.matrix(dist(coords_matrix))

# Parâmetros do problema p-centro
p <-7  # Número de centros a serem selecionados
n <- nrow(bairros_itacoatiara) # Número de bairros (demandas e potenciais locais para centros)

# Implementar o modelo p-centro

# Variáveis de decisão: y_ij (n*n), x_j (n), w (1)
# Função objetivo: minimizar w (distância máxima de cobertura)
obj <- c(rep(0, n * n), rep(0, n), 1)

# Determinar o número total de restrições
num_constraints <- n + 1 + n * n + n * n

# Inicializar matriz A, vetor b e vetor dir (direções das restrições)
A <- matrix(0, nrow = num_constraints, ncol = n * n + n + 1)
b <- numeric(num_constraints)
dir <- character(num_constraints)

current_row <- 0 # Contador para a linha atual na matriz A

# Restrição 1 (A1): Cada bairro i deve ser atendido por exatamente um centro j
for (i in 1:n) {
  current_row <- current_row + 1
  A[current_row, ((i - 1) * n + 1):(i * n)] <- 1 # Coeficientes para y_ij
  b[current_row] <- 1
  dir[current_row] <- "="
}

# Restrição 2 (A2): Exatamente p centros devem ser selecionados
for (i in 1:n) {
  current_row <- current_row + 1
  A[current_row, (n * n + 1):(n * n + n)] <- 1 # Coeficientes para x_j
  b[current_row] <- p
  dir[current_row] <- "="
}

# Restrição 3 (A3): Y_ij <= X_j 
for (i in 1:n) {
  for (j in 1:n) {
    current_row <- current_row + 1
    idx_y_ij <- (i - 1) * n + j 
    idx_x_j <- n * n + j       
    A[current_row, idx_y_ij] <- 1
    A[current_row, idx_x_j] <- -1
    b[current_row] <- 0
    dir[current_row] <- "<="
  }
}

# Restrição 4 (A4): d_ij * Y_ij <= w 
for (i in 1:n) {
  for (j in 1:n) {
    current_row <- current_row + 1
    idx_y_ij <- (i - 1) * n + j       
    idx_w <- n * n + n + 1            
    A[current_row, idx_y_ij] <- dist_matrix[i, j] 
    A[current_row, idx_w] <- -1                   
    b[current_row] <- 0
    dir[current_row] <- "<="
  }
}

# Resolver o problema de otimização
# indices_bin são as variáveis y_ij e x_j
indices_bin <- 1:(n * n + n) 

resultado_p_centro <- lp("min", obj, A, dir, b, binary = indices_bin)


# Processar e mapear os resultados do p-centro 
bairros_itacoatiara <- bairros_itacoatiara %>%
  mutate(p_centro_selected = "Não Selecionado")

if (resultado_p_centro$status == 0) {
  cat("Solução p-centro encontrada!\n")
  solucao_p_centro <- resultado_p_centro$solution
  
  # Variável renomeada: w_otimo para distancia_max_cobertura_otima
  distancia_max_cobertura_otima <- solucao_p_centro[n * n + n + 1] 
  
  cat("Distância máxima de cobertura ótima:", distancia_max_cobertura_otima, "\n")
  
  # As variáveis x_j (seleção de centro) começam em n*n + 1 e vão até n*n + n
  postos_otimos_indices_centro <- which(solucao_p_centro[(n * n + 1):(n * n + n)] == 1)
  bairros_otimos_centro <- bairros_itacoatiara$name_neighborhood[postos_otimos_indices_centro]
  cat("Os locais ótimos (p-centro) são nos bairros:\n", paste(bairros_otimos_centro, collapse = ", "), "\n")
  bairros_itacoatiara$p_centro_selected[postos_otimos_indices_centro] <- "Selecionado (Ótimo)"
} else {
  cat("Não foi possível encontrar uma solução ótima para o p-centro. Status:", resultado_p_centro$status, "\n")
}

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
    subtitle = "Bairros preenchidos pela categoria"
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