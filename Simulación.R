rm(list=ls())

suppressMessages(suppressWarnings({
  packages <- c("sand", "dplyr", "jsonlite", "tidyr", "igraph", "corrplot", "gplots", "xtable","viridis")
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}))

library(igraph)



#-------------------------------------------------------------------------------
#Obtener grafos con diferentes características estructurales

# Función para calcular las métricas solicitadas y devolver un data.frame
g_metr <- function(graph, name) {
  densidad <- graph.density(graph)
  transitividad <- transitivity(graph, type = "global")
  grado_promedio <- mean(degree(graph))
  sd_grado <- sd(degree(graph))
  
  if (is_connected(graph)) {
    distancia_media <- mean_distance(graph)
  } else {
    distancia_media <- mean_distance(graph, unconnected = TRUE)
  }
  
  # Devuelve un data.frame con una fila de métricas
  return(data.frame(
    Grafo = name,
    Densidad = densidad,
    Transitividad = transitividad,
    GradoMedio = grado_promedio,
    SDGrado = sd_grado,
    DistanciaMedia = distancia_media
  ))
}


#Grafos
# Estructura 2: Grafo Small-World
set.seed(784)
g_smallworld <- sample_smallworld(dim = 1, size = 24, nei = 3, p = 0.1)
plot(g_smallworld, layout = layout_with_fr, vertex.size = 5, vertex.label = NA, main = "Grafo Small-World")


# Estructura 3: Grafo Barabási-Albert
set.seed(524)
g_barabasi <- barabasi.game(n = 24, m = 2, directed = FALSE)
plot(g_barabasi, layout = layout_with_fr, vertex.size = 10, main = "Grafo Barabási-Albert")


# Estructura 4: Árbol Barabási-Albert
set.seed(876)
g_tree_barabasi <- sample_pa(n = 24, m = 1, directed = FALSE)
plot(g_tree_barabasi, layout = layout_with_fr, vertex.size = 5, vertex.label = NA, main = "Árbol Barabási-Albert")

# Estructura 5: Grafo Aleatorio (Erdos-Rényi)
set.seed(75)
g_random <- sample_gnp(n = 24, p = 0.1)
plot(g_random, layout = layout_with_fr, vertex.size = 5, vertex.label = NA, main = "Grafo Erdos-Rényi")

# Estructura 6: Grafo Strike
data(strike)
upgrade_graph(strike)
plot(strike, vertex.size = 5, vertex.label = NA, main = "Grafo Strike")


#-------------------------------------------------------------------------------
# Caracterización de la estructura
g_list <- list(MundoPequeño = g_smallworld, BarabásiAlbert = g_barabasi, Huelga = strike, GrafoAleatorio = g_random, ÁrbolBarabásiAlbert = g_tree_barabasi)
results <- lapply(names(g_list), function(name) g_metr(g_list[[name]], name))
output_matrix <- do.call(rbind, results)
output_matrix[, -1] <- round(output_matrix[, -1], 3)
print(output_matrix)

#-------------------------------------------------------------------------------
#Simulación

# Función para el análisis y simulación en base a un grafo
ACE_MCsampling <- function(graph, nv_values, m_values, base_repetitions, O_values) {
  results <- list()
  
  # Iterar sobre los valores de nv y m
  for (nv in nv_values) {
    for (m in m_values) {
      # Ajustar las repeticiones proporcionalmente al valor de m
      if (m == min(m_values)){
        repetitions <- base_repetitions
      } else{
        repetitions <- base_repetitions *m   
      }
      
      A1 <- as_adjacency_matrix(graph)
      for (n in repetitions) {
        set.seed(41)
        I11 <- matrix(, nrow = nv, ncol = n)
        I10 <- matrix(, nrow = nv, ncol = n)
        I01 <- matrix(, nrow = nv, ncol = n)
        I00 <- matrix(, nrow = nv, ncol = n)
        
        for (i in 1:n) {
          z <- rep(0, nv)
          sample <- sample(1:nv, m, replace = FALSE)
          z[sample] <- 1
          sample_nbrs <- as.numeric(z %*% A1 > 0)
          I11[, i] <- z * sample_nbrs
          I10[, i] <- z * (1 - sample_nbrs)
          I01[, i] <- (1 - z) * sample_nbrs
          I00[, i] <- (1 - z) * (1 - sample_nbrs)
        }
        
        I11.11 <- I11 %*% t(I11) / n
        I10.10 <- I10 %*% t(I10) / n
        I01.01 <- I01 %*% t(I01) / n
        I00.00 <- I00 %*% t(I00) / n
        
        Obar.c11 <- numeric()
        Obar.c10 <- numeric()
        Obar.c01 <- numeric()
        Obar.c00 <- numeric()
        
        c11_values <- numeric()
        c10_values <- numeric()
        c01_values <- numeric()
        c00_values <- numeric()
        
        for (i in 1:n) {
          z <- rep(0, nv)
          sample <- sample(1:nv, m, replace = FALSE)
          z[sample] <- 1
          sample_nbrs <- as.numeric(z %*% A1 > 0)
          
          c11 <- z * sample_nbrs
          c10 <- z * (1 - sample_nbrs)
          c01 <- (1 - z) * sample_nbrs
          c00 <- (1 - z) * (1 - sample_nbrs)
          
          Obar.c11 <- c(Obar.c11, O_values[1] * mean(c11 / diag(I11.11)))
          Obar.c10 <- c(Obar.c10, O_values[2] * mean(c10 / diag(I10.10)))
          Obar.c01 <- c(Obar.c01, O_values[3] * mean(c01 / diag(I01.01)))
          Obar.c00 <- c(Obar.c00, O_values[4] * mean(c00 / diag(I00.00)))
          
          # Guardar los valores de c11, c10, c01, y c00 para cada iteración
          c11_values <- c(c11_values, sum(c11))
          c10_values <- c(c10_values, sum(c10))
          c01_values <- c(c01_values, sum(c01))
          c00_values <- c(c00_values, sum(c00))
        }
        
        Obar <- list(Obar.c11, Obar.c10, Obar.c01, Obar.c00)
        ACE <- list(Obar.c11-Obar.c00, Obar.c10-Obar.c00, Obar.c01-Obar.c00)
        c_values <- list(c11_values, c10_values, c01_values, c00_values)
        
        means <- sapply(ACE, mean)
        sds <- sapply(ACE, sd)
        cvs <- sapply(ACE, sd) / sapply(ACE, mean)
        
        results[[paste("nv=", nv, "m=", m, "reps=", n)]] <- list(
          means = means, sds = sds, cvs = cvs, c_values = c_values, ACE = ACE, Obar = Obar, I11=I11.11, I10=I10.10, I01=I01.01, I00=I00.00
        )
      }
    }
  }
  
  return(results)
}



#1. Considerar múltiples repeticiones- número fijo de tratados
# Parámetros comunes
nv_values <- c(24)
m_values <- c(4, 8, 12) # Diferentes números de tratados
O_values <- c(10, 7, 5, 1) # Pagos
base_repetitions <- c(1000, 5000, 10000) # Repeticiones base para m = 4


g_list <- list(
  "MundoPequeño" = g_smallworld,
  "BarabásiAlbert" = g_barabasi,  
  "Huelga" = strike,
  "GrafoAleatorio" = g_random,
  "ÁrbolBarabásiAlbert" = g_tree_barabasi
)

# Calcular el efecto causal
results <- list()

for (name in names(g_list)) {
  cat("Grafo:", name, "\n")
  graph <- g_list[[name]]
  
  results[[name]] <- ACE_MCsampling(graph, nv_values, m_values, repetitions, O_values)
  result <- results[[name]][[as.character(max(repetitions))]]
  
}

#-------------------------------------------------------------------------------
#Resultados
table_results <- list()

# Efecto Global: Obar.c11 - Obar.c00 
for (i in 1:length(g_list)) {
  grafo_name <- names(g_list)[i]
  media <- c(
    results[[i]][[1]]$means[1], results[[i]][[2]]$means[1], results[[i]][[3]]$means[1], results[[i]][[4]]$means[1], 
    results[[i]][[5]]$means[1], results[[i]][[6]]$means[1], results[[i]][[7]]$means[1], results[[i]][[8]]$means[1], 
    results[[i]][[9]]$means[1])
  se <- c(
    results[[i]][[1]]$sds[1], results[[i]][[2]]$sds[1], results[[i]][[3]]$sds[1], results[[i]][[4]]$sds[1],
    results[[i]][[5]]$sds[1], results[[i]][[6]]$sds[1], results[[i]][[7]]$sds[1], results[[i]][[8]]$sds[1], 
    results[[i]][[9]]$sds[1])
  cv <- c(
    results[[i]][[1]]$cvs[1], results[[i]][[2]]$cvs[1], results[[i]][[3]]$cvs[1], results[[i]][[4]]$cvs[1], 
    results[[i]][[5]]$cvs[1], results[[i]][[6]]$cvs[1], results[[i]][[7]]$cvs[1], results[[i]][[8]]$cvs[1], 
    results[[i]][[9]]$cvs[1])
  
  # Almacenar los resultados en table_results
  table_results[[grafo_name]] <- rbind(media, se, cv)
}
output_matrix_ace <- do.call(rbind, table_results)
print(round(output_matrix_ace,2))

# Efecto directo: Obar.c10 - Obar.c00 
for (i in 1:length(g_list)) {
  grafo_name <- names(g_list)[i]
  media <- c(
    results[[i]][[1]]$means[2], results[[i]][[2]]$means[2], results[[i]][[3]]$means[2], results[[i]][[4]]$means[2], 
    results[[i]][[5]]$means[2], results[[i]][[6]]$means[2], results[[i]][[7]]$means[2], results[[i]][[8]]$means[2], 
    results[[i]][[9]]$means[2])
  se <- c(
    results[[i]][[1]]$sds[2], results[[i]][[2]]$sds[2], results[[i]][[3]]$sds[2], results[[i]][[4]]$sds[2],
    results[[i]][[5]]$sds[2], results[[i]][[6]]$sds[2], results[[i]][[7]]$sds[2], results[[i]][[8]]$sds[2], 
    results[[i]][[9]]$sds[2])
  cv <- c(
    results[[i]][[1]]$cvs[2], results[[i]][[2]]$cvs[2], results[[i]][[3]]$cvs[2], results[[i]][[4]]$cvs[2], 
    results[[i]][[5]]$cvs[2], results[[i]][[6]]$cvs[2], results[[i]][[7]]$cvs[2], results[[i]][[8]]$cvs[2], 
    results[[i]][[9]]$cvs[2])
  
  # Almacenar los resultados en table_results
  table_results[[grafo_name]] <- rbind(media, se, cv)
}
output_matrix_ace <- do.call(rbind, table_results)
print(round(output_matrix_ace,2))

# Efecto indirectos: Obar.c01 - Obar.c00 
for (i in 1:length(g_list)) {
  grafo_name <- names(g_list)[i]
  media <- c(
    results[[i]][[1]]$means[3], results[[i]][[2]]$means[3], results[[i]][[3]]$means[3], results[[i]][[4]]$means[3], 
    results[[i]][[5]]$means[3], results[[i]][[6]]$means[3], results[[i]][[7]]$means[3], results[[i]][[8]]$means[3], 
    results[[i]][[9]]$means[3])
  se <- c(
    results[[i]][[1]]$sds[3], results[[i]][[2]]$sds[3], results[[i]][[3]]$sds[3], results[[i]][[4]]$sds[3],
    results[[i]][[5]]$sds[3], results[[i]][[6]]$sds[3], results[[i]][[7]]$sds[3], results[[i]][[8]]$sds[3], 
    results[[i]][[9]]$sds[3])
  cv <- c(
    results[[i]][[1]]$cvs[3], results[[i]][[2]]$cvs[3], results[[i]][[3]]$cvs[3], results[[i]][[4]]$cvs[3], 
    results[[i]][[5]]$cvs[3], results[[i]][[6]]$cvs[3], results[[i]][[7]]$cvs[3], results[[i]][[8]]$cvs[3], 
    results[[i]][[9]]$cvs[3])
  
  # Almacenar los resultados en table_results
  table_results[[grafo_name]] <- rbind(media, se, cv)
}
output_matrix_ace <- do.call(rbind, table_results)
print(round(output_matrix_ace,2))


#-------------------------------------------------------------------------------
# Distribución c11, c10, c01, c00- Barabási-Albert, 120.000, T==12
boxplot(results[[2]][[9]]$c_values[[1]], results[[2]][[9]]$c_values[[2]], results[[2]][[9]]$c_values[[3]], results[[2]][[9]]$c_values[[4]], 
        main = "Nodos según nival de exposición\nGrafo: Barabási-Albert; T: 12; Iteraciones: 120,000",
        names = c("c11", "c10", "c01", "c00"))
table(diag(round(results[[2]][[9]]$I00, 2)))
table(degree(g_barabasi))

boxplot(results[[1]][[3]]$ACE[[1]], results[[1]][[3]]$ACE[[2]], results[[1]][[3]]$ACE[[3]], 
        main = "Nodos según nival de exposición\nGrafo: Barabási-Albert; T: 12; Iteraciones: 120,000",
        names = c("Obar.c11- Obar.c00", "Obar.c10- Obar.c00", "Obar.c01- Obar.c00"))
boxplot(results[[4]][[3]]$ACE[[1]], results[[4]][[3]]$ACE[[2]], results[[4]][[3]]$ACE[[2]], 
        main = "Nodos según nival de exposición\nGrafo: Grafo Aleatorio; T: 12; Iteraciones: 120,000",
        names = c("Obar.c11- Obar.c00", "Obar.c10- Obar.c00", "Obar.c01- Obar.c00"))
boxplot(results[[5]][[3]]$ACE[[1]], results[[5]][[3]]$ACE[[2]], results[[5]][[3]]$ACE[[2]], 
        main  = "Nodos según nival de exposición\nGrafo: Árbol Barabási-Albert; T: 12; Iteraciones: 120,000",
        names = c("Obar.c11- Obar.c00", "Obar.c10- Obar.c00", "Obar.c01- Obar.c00"))


