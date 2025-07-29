library(lpSolve)
library(geobr)
library(sf)
library(ggplot2)
library(dplyr)
library(units) # Necessário para drop_units

# Definir o raio máximo de cobertura (em KM)
R_max_km <- 2 # Exemplo: 2 quilômetros. Ajuste conforme necessário.
cat("Raio máximo de cobertura (R_max_km):", R_max_km, "km\n")

# Dados dos bairros de Itacoatiara
bairros <- read_neighborhood(year = 2022)
bairros_itacoatiara <- bairros %>%
  filter(code_muni == 1301902) %>%
  select(code_neighborhood, name_neighborhood, geom) %>%
  mutate(label_id = row_number())

if (nrow(bairros_itacoatiara) == 0) {
  stop("Não foram encontrados bairros para Itacoatiara-AM")
}

# Calcular a matriz de distâncias entre os centroides dos bairros
centroids <- st_centroid(bairros_itacoatiara$geom)
dist_matrix <- st_distance(centroids, centroids) 
dist_matrix_km <- drop_units(dist_matrix) / 1000 

# Implementar o modelo de Cobertura de Conjuntos
n <- nrow(bairros_itacoatiara)

# Variáveis de decisão (x_j): 1 se um posto é instalado no bairro j, 0 caso contrário
# Função objetivo: minimizar o número de postos (sum(x_j))
obj <- rep(1, n)

# Matriz de Cobertura (A): A_ij = 1 se o bairro i pode ser coberto por um posto em j
# (ou seja, dist(i,j) <= R_max_km)
matrix_cobertura <- ifelse(dist_matrix_km <= R_max_km, 1, 0)

# Restrições: Cada bairro deve ser coberto por pelo menos um posto
A <- matrix_cobertura 
b <- rep(1, n)        
dir <- rep(">=", n)   

# Resolver o problema de otimização
cat("Tentando resolver o problema de Cobertura de Conjuntos com lpSolve...\n")
resultado_cobertura <- lp("min", obj, A, dir, b, all.bin = TRUE) 
cat("lpSolve finalizado. Verificando o status...\n")

# Processar e mapear os resultados da Cobertura de Conjuntos 
bairros_itacoatiara <- bairros_itacoatiara %>%
  mutate(cobertura_selected = "Não Selecionado") 

if (resultado_cobertura$status == 0) {
  cat("Solução de Cobertura de Conjuntos encontrada!\n")
  solucao_cobertura <- resultado_cobertura$solution
  postos_otimos_indices_cobertura <- which(solucao_cobertura == 1)
  bairros_otimos_cobertura_nomes <- bairros_itacoatiara$name_neighborhood[postos_otimos_indices_cobertura]
  cat("Os locais ótimos (Cobertura de Conjuntos) são nos bairros:\n", paste(bairros_otimos_cobertura_nomes, collapse = ", "), "\n")
  bairros_itacoatiara$cobertura_selected[postos_otimos_indices_cobertura] <- "Selecionado (Ótimo)" # Ajuste para "Ótimo"
} else {
  cat("Não foi possível encontrar uma solução ótima para Cobertura de Conjuntos. Status do solver:", resultado_cobertura$status, "\n")
  cat("Possíveis razões para o status (veja a documentação do lpSolve para detalhes):\n")
  cat(" - 2: O problema é inviável (i.e., com o R_max_km dado, alguns bairros não podem ser cobertos).\n")
}

#Carregar e filtrar dados das UBS existentes
servicos_saude_brasil <- read_health_facilities() %>%
  st_as_sf()
ubs_nomes_selecionados <- c(
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
    no_fantasia %in% ubs_nomes_selecionados
  ) %>%
  rename(name = no_fantasia) %>%
  select(name, no_bairro, geom)

if (nrow(servicos_saude) == 0) {
  warning("Nenhuma UBS encontrada com os referidos nomes em Itacoatiara-AM. O mapa pode não mostrar 'Existente' ou 'Ótimo e Existente'.")
} else {
  cat("\nUBS existentes selecionadas e consideradas:\n")
  print(servicos_saude %>% st_drop_geometry() %>% select(name, no_bairro))
}


# Realizar a junção espacial para associar serviços de saúde aos bairros
servicos_saude_bairros_spatial <- st_join(servicos_saude, bairros_itacoatiara, join = st_within, left = FALSE) %>%
  select(name, label_id)

bairros_com_servicos <- servicos_saude_bairros_spatial %>%
  distinct(label_id) %>%
  pull(label_id)

# Classificar os bairros para mapeamento 
bairros_itacoatiara <- bairros_itacoatiara %>%
  mutate(
    tem_posto_existente = ifelse(label_id %in% bairros_com_servicos, "Existente", "Não Existente"),
    category = case_when( # Nome da coluna para consistência
      cobertura_selected == "Selecionado (Ótimo)" & tem_posto_existente == "Existente" ~ "Atual/Ótima",
      cobertura_selected == "Selecionado (Ótimo)" & tem_posto_existente == "Não Existente" ~ "Localização Ótima",
      tem_posto_existente == "Existente" & cobertura_selected == "Não Selecionado" ~ "Localização UBS Atual",
      TRUE ~ "Nenhum"
    )
  )

# Definir cores para o mapa com as novas legendas
cores_mapa <- c( # Nome do vetor para consistência
  "Atual/Ótima" = "purple",
  "Localização Ótima" = "blue",
  "Localização UBS Atual" = "green",
  "Nenhum" = "lightgrey" # Cor para bairros sem nenhum status
)

# Gerar o mapa 
mapa_bairros_numerados_colorido <- ggplot() + # Nome do objeto para consistência
  geom_sf(data = bairros_itacoatiara, aes(fill = category), color = "black", linewidth = 0.2) +
  # >>> MODIFICAÇÃO AQUI: Removido o filtro para mostrar todos os labels <<<
  geom_sf_text(data = bairros_itacoatiara, # Apenas especificar os dados aqui, sem filtro
               aes(label = label_id, geometry = st_centroid(geom)),
               size = 2, color = "black", check_overlap = TRUE, fontface = "bold") +
  scale_fill_manual(values = cores_mapa, name = "Legenda") +
  labs(
    title = paste0("Localização de Pontos Ótimos (Cobertura - ", R_max_km, "km) vs. UBS Selecionadas em Itacoatiara"),
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