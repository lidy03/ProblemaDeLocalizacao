library(lpSolve)
library(geobr)
library(sf)
library(ggplot2)
library(dplyr)
library(units) 

# Raio Máximo de Cobertura (km) 
R_max_km <- 2 
cat("Raio máximo de cobertura (R_max_km):", R_max_km, "km\n")

# Dados dos bairros de Itacoatiara 
bairros <- read_neighborhood(year = 2022)
bairros_itacoatiara <- bairros %>%
  filter(code_muni == 1301902) %>%
  select(code_neighborhood, name_neighborhood, geom) %>%
  mutate(label_id = row_number()) 
n <- nrow(bairros_itacoatiara)

if (n == 0) {
  stop("Não foram encontrados bairros para Itacoatiara-AM")
}

# Calcular a matriz de distâncias entre os centroides dos bairros
centroids <- st_centroid(bairros_itacoatiara$geom)
dist_matrix <- st_distance(centroids, centroids) 
dist_matrix_km <- drop_units(dist_matrix) / 1000 

# Implementar o modelo 

obj <- rep(1, n)

matrix_cobertura <- ifelse(dist_matrix_km <= R_max_km, 1, 0)

# Restrições: Cada bairro deve ser coberto por pelo menos um posto
A <- matrix_cobertura
b <- rep(1, n)        
dir <- rep(">=", n)   

# Resolver o problema de otimização 
cat("Tentando resolver o problema de Cobertura de Conjuntos...\n")
resultado_cobertura <- lp("min", obj, A, dir, b, all.bin = TRUE) 
cat("lpSolve finalizado. Verificando a solução...\n")

# Processar os resultados da Cobertura de Conjuntos
bairros_itacoatiara <- bairros_itacoatiara %>%
  mutate(cobertura_selected = "Não Selecionado") 

if (resultado_cobertura$status == 0) {
  solucao_cobertura <- resultado_cobertura$solution
  postos_otimos_indices_cobertura <- which(solucao_cobertura == 1)
  
  cat("\nSolução de Cobertura de Conjuntos encontrada\n")
  cat("Número mínimo de postos necessários:", resultado_cobertura$objval, "\n")
  
  bairros_otimos_cobertura_nomes <- bairros_itacoatiara$name_neighborhood[postos_otimos_indices_cobertura]
  cat("Os locais ótimos (Cobertura de Conjuntos) são nos bairros:\n", paste(bairros_otimos_cobertura_nomes, collapse = ", "), "\n")
  
  bairros_itacoatiara$cobertura_selected[postos_otimos_indices_cobertura] <- "Selecionado (Ótimo)"
} else {
  cat("Não foi possível encontrar uma solução ótima para Cobertura de Conjuntos. Status:", resultado_cobertura$status, "\n")
}

# Mapa para visualização dos pontos ótimos 
mapa_otimizacao <- ggplot() +
  geom_sf(data = bairros_itacoatiara, 
          aes(fill = cobertura_selected), 
          color = "black", 
          linewidth = 0.2) +
  geom_sf_text(data = bairros_itacoatiara,
               aes(label = label_id, geometry = st_centroid(geom)),
               size = 2, color = "black", check_overlap = TRUE, fontface = "bold") +
  
  scale_fill_manual(values = c("Selecionado (Ótimo)" = "purple", "Não Selecionado" = "lightgrey"),
                    name = "Status de Seleção") +
  labs(
    title = paste0("Localização Ótima (Cobertura de Conjuntos - ", R_max_km, "km)"),
    subtitle = "Bairros de Itacoatiara-AM"
  ) +
  theme_minimal()

print(mapa_otimizacao)

#Imprimir a legenda numérica dos bairros no console
legenda_numerica_bairros <- bairros_itacoatiara %>%
  st_drop_geometry() %>%
  select(ID = label_id, Bairro = name_neighborhood, Cod_Bairro_IBGE = code_neighborhood) %>%
  arrange(ID)

print("\n--- Legenda dos Bairros ---")
print(legenda_numerica_bairros)
cat("--------------------------\n")