library(geobr)
library(sf)
library(ggplot2)
library(dplyr)
library(lpSolve)

#Dados dos bairros de Itacoatiara
bairros <- read_neighborhood(year = 2022)

bairros_itacoatiara <- bairros %>%
  filter(code_muni == 1301902) %>%
  select(code_neighborhood, name_neighborhood, geom) %>%
  mutate(label_id = row_number()) 

if (nrow(bairros_itacoatiara) == 0) {
  stop("Não foram encontrados bairros para Itacoatiara-AM")
}

# Matriz de distâncias entre os centroides dos bairros
centroides <- st_centroid(bairros_itacoatiara$geom)
dist_matrix <- as.matrix(dist(st_coordinates(centroides)))

n <- nrow(bairros_itacoatiara) 
p <- 7                      

cat("Resolvendo o problema p-mediana com p =", p, "...\n")

#Implementar o modelo

# Variáveis de decisão: y_ij (n*n) + x_j (n)
# Função objetivo: minimizar a soma das distâncias (dist_matrix) ponderadas por y_ij (obj)
# Coeficientes para y_ij (dist_matrix) e coeficientes zero para x_j
obj <- c(as.vector(dist_matrix), rep(0, n))

# Restrições A1: Cada bairro i deve ser atendido por exatamente um centro j (sum(y_ij para j) = 1)
# Matriz A1 (n linhas x (n*n + n) colunas)
A1 <- matrix(0, n, n * n + n)
for (i in 1:n) {
  A1[i, ((i - 1) * n + 1):(i * n)] <- 1
}
b1 <- rep(1, n)

# Restrições A2: Se o bairro i é atendido pelo bairro j (y_ij=1), o bairro j deve ter um centro (x_j=1). (y_ij <= x_j)
# Matriz A2 (n*n linhas x (n*n + n) colunas)
A2 <- matrix(0, n * n, n * n + n)
for (i in 1:n) {
  for (j in 1:n) {
    idx <- (i - 1) * n + j 
    A2[idx, idx] <- 1
    A2[idx, n * n + j] <- -1 
  }
}
b2 <- rep(0, n * n)

# Restrição A3: Exatamente p centros devem ser selecionados (sum(x_j) = p)
# Matriz A3 (1 linha x (n*n + n) colunas)
A3 <- matrix(0, 1, n * n + n)
A3[1, (n * n + 1):(n * n + n)] <- 1 
b3 <- p

# Montar a matriz de restrições A, vetor de direções e vetor b
A <- rbind(A1, A2, A3)
dir <- c(rep("=", n), rep("<=", n * n), "=")
b <- c(b1, b2, b3)

# Resolver o problema de otimização (todas as variáveis y_ij e x_j são binárias)
resultado_p_mediana <- lp("min", obj, A, dir, b, all.bin = TRUE)

# Processar os resultados

# Inicializar coluna de seleção ótima
bairros_itacoatiara <- bairros_itacoatiara %>%
  mutate(p_mediana_selected = "Não Selecionado")

if (resultado_p_mediana$status == 0) {
  cat("\nSolução p-mediana encontrada!\n")
  
  solucao_p_mediana_x <- resultado_p_mediana$solution[(n * n + 1):(n * n + n)]
  
  postos_otimos_indices <- which(solucao_p_mediana_x == 1)
  bairros_otimos <- bairros_itacoatiara$name_neighborhood[postos_otimos_indices]
  
  cat("A soma mínima das distâncias de cobertura é:", resultado_p_mediana$objval, "\n")
  cat("Os locais ótimos (p-mediana) são nos bairros:\n", paste(bairros_otimos, collapse = ", "), "\n")
  
  bairros_itacoatiara$p_mediana_selected[postos_otimos_indices] <- "Selecionado (Ótimo)"
} else {
  cat("Não foi possível encontrar uma solução ótima para o p-mediana. Status:", resultado_p_mediana$status, "\n")
}

# Mapa simples para visualização dos pontos ótimos
mapa_otimizacao_pmediana <- ggplot() +
  geom_sf(data = bairros_itacoatiara, 
          aes(fill = p_mediana_selected), 
          color = "black", 
          linewidth = 0.2) +
  geom_sf_text(data = bairros_itacoatiara,
               aes(label = label_id, geometry = st_centroid(geom)),
               size = 2, color = "black", check_overlap = TRUE, fontface = "bold") +
  
  scale_fill_manual(values = c("Selecionado (Ótimo)" = "purple", "Não Selecionado" = "lightgrey"),
                    name = "Status de Seleção") +
  labs(
    title = paste0("Localização Ótima (p-Mediana, p=", p, ")"),
    subtitle = "Bairros de Itacoatiara-AM"
  ) +
  theme_minimal()

print(mapa_otimizacao_pmediana)

# Imprimir a legenda numérica dos bairros no console
legenda_numerica <- bairros_itacoatiara %>%
  st_drop_geometry() %>%
  select(ID = label_id, Bairro = name_neighborhood) %>%
  arrange(ID)

print("\n--- Legenda dos Bairros ---")
print(legenda_numerica)
cat("--------------------------\n")