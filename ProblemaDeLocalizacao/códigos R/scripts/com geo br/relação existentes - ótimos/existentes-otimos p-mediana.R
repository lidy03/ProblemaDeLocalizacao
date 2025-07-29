library(geobr)
library(sf)
library(ggplot2)
library(dplyr)
library(lpSolve)

# Baixar os dados dos bairros do Brasil
bairros <- read_neighborhood(year = 2022)

bairros_itacoatiara <- bairros %>%
  filter(code_muni == 1301902) %>%
  select(code_neighborhood, name_neighborhood, geom) %>%
  mutate(label_id = row_number()) 

if (nrow(bairros_itacoatiara) == 0) {
  stop("Não foram encontrados bairros para Itacoatiara-AM")
}

# Calcular a matriz de distâncias entre os centroides dos bairros
centroides <- st_centroid(bairros_itacoatiara$geom)
dist_matrix <- as.matrix(dist(st_coordinates(centroides)))
n <- nrow(bairros_itacoatiara)
p <- 7 

# Implementar o modelo
obj <- c(as.vector(dist_matrix), rep(0, n))
A1 <- matrix(0, n, n * n + n)
for (i in 1:n) {
  A1[i, ((i - 1) * n + 1):(i * n)] <- 1
}
b1 <- rep(1, n)
A2 <- matrix(0, n * n, n * n + n)
for (i in 1:n) {
  for (j in 1:n) {
    idx <- (i - 1) * n + j
    A2[idx, idx] <- 1
    A2[idx, n * n + j] <- -1
  }
}
b2 <- rep(0, n * n)
A3 <- matrix(0, 1, n * n + n)
A3[1, (n * n + 1):(n * n + n)] <- 1
b3 <- p
A <- rbind(A1, A2, A3)
dir <- c(rep("=", n), rep("<=", n * n), "=")
b <- c(b1, b2, b3)
resultado_p_mediana <- lp("min", obj, A, dir, b, all.bin = TRUE)

# Inicializar coluna de seleção ótima
bairros_itacoatiara <- bairros_itacoatiara %>%
  mutate(p_mediana_selected = "Não Selecionado")

if (resultado_p_mediana$status == 0) {
  cat("Solução p-mediana encontrada!\n")
  solucao_p_mediana <- resultado_p_mediana$solution[(n * n + 1):(n * n + n)]
  postos_otimos_indices <- which(solucao_p_mediana == 1)
  bairros_otimos <- bairros_itacoatiara$name_neighborhood[postos_otimos_indices]
  cat("Os locais ótimos (p-mediana) são nos bairros:\n", paste(bairros_otimos, collapse = ", "), "\n")
  bairros_itacoatiara$p_mediana_selected[postos_otimos_indices] <- "Selecionado (Ótimo)"
} else {
  cat("Não foi possível encontrar uma solução ótima para o p-mediana.\n")
}

servicos_saude_brasil <- read_health_facilities() %>%
  st_as_sf()

servicos_saude <- servicos_saude_brasil %>%
  filter(
    abbrev_state == "AM",
    code_muni == 1301902,
    no_fantasia %in% ubs_nomes
  ) %>%
  rename(name = no_fantasia) %>%
  select(name, geom)

# Realizar a junção espacial para associar serviços de saúde aos bairros
servicos_saude_bairros_spatial <- st_join(servicos_saude, bairros_itacoatiara, join = st_within, left = FALSE) %>%
  select(name, label_id) # Usar label_id para vincular

bairros_com_servicos <- servicos_saude_bairros_spatial %>%
  distinct(label_id) %>%
  pull(label_id)

bairros_itacoatiara <- bairros_itacoatiara %>%
  mutate(
    tem_posto_existente = ifelse(label_id %in% bairros_com_servicos, "Existente", "Não Existente"),
    category = case_when(
      p_mediana_selected == "Selecionado (Ótimo)" & tem_posto_existente == "Existente" ~ "Atual/Ótimo",
      p_mediana_selected == "Selecionado (Ótimo)" & tem_posto_existente == "Não Existente" ~ "Localização ótima",
      tem_posto_existente == "Existente" & p_mediana_selected == "Não Selecionado" ~ "Localização UBS existente",
      TRUE ~ "Nenhum"
    )
  )

cores_mapa <- c(
  "Ótimo e Existente" = "purple",
  "Apenas Ótimo" = "blue",
  "Apenas Existente" = "green",
  "Nenhum" = "white"
)

mapa_bairros_numerados_colorido <- ggplot() +
  geom_sf(data = bairros_itacoatiara, aes(fill = category), color = "black") +
  geom_sf_text(data = bairros_itacoatiara, aes(label = label_id), size = 2, color = "black") + 
  scale_fill_manual(values = cores_mapa, name = "Legenda") + 
  labs(
    title = "Localização de Pontos Ótimos (p-mediana) vs. Serviços de Saúde em Itacoatiara",
    subtitle = "Bairros preenchidos pela categoria"
  ) +
  theme_minimal()

print(mapa_bairros_numerados_colorido)

legenda_numerica <- bairros_itacoatiara %>%
  select(ID = label_id, Bairro = name_neighborhood) %>%
  arrange(ID)

print("\nLegenda dos Bairros:")
print(legenda_numerica)
