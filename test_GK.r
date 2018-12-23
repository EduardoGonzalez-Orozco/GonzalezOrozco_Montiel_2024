
### Código con librerías para probar el estadìstico tau de Goodman-Kruskal


#Librerías:
library(DescTools)
library(genetics)

# Input genotipos 
g1 <-  genotype( c("T/T","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A","T/T","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A"))
g2 <-  genotype( c("G/C","G/G","G/G","G/G","G/G","G/G","G/G","G/G","G/G","G/C","G/G","G/G","G/G","G/G","G/G","G/G","G/G","G/G","G/G","G/G"))

#Estimación utilizando estadísticos clásicos 
lin <-LD(g1,g2)

##Output de resultados
lin$`R^2`
lin$`D'`
lin$D


##Input para tau de Goodman Kruskal
gg1 <-  c("T/T","A/A","A/A","T/T","A/A","A/A","A/A","A/A","T/T","T/T","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A","A/A")
gg2 <-   c("G/C","G/G","G/G","G/C","G/G","G/G","G/G","G/G","G/C","G/C","G/G","G/G","G/G","G/G","G/G","G/G","G/G","G/G","G/G","G/G")

# Estimaación utilizando Goodman_kruskal tau
GoodmanKruskalTau(gg1,gg2,direction = "column")

## Calculo y output de resultados
GoodmanKruskalTau(gg1,gg2)




