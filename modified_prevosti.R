## Die Funktion ModProv berechnet eine genetische Distanzmatrix zwischen Individuen. Sie benoetigt als 
## Parameter einen Dataframe oder eine Matrix mit Genotypen in Zeilen und Markern in Spalten, mit jeweils
## zwei Spalten pro Marker (eine pro Allel) in denen das Vorhandensein des jeweiligen Allels mit 1 
## angegeben wird, bzw die Abwesenheit des Allels mit 0. Möglich sind also pro Marker und Genotyp die 
## Kombinationen 10, 01 und 11. Per Definition haben zwei gleiche Homozygote zueinander die Distanz 0, 
## zwei unterschiedliche Homozygote haben die Distanz 1 und Heterozygote zu Homozygoten, sowie zwei Hetero-
## zygote zueinander haben die Distanz 0.5.
##########################################################################################################

## ModProv v2.0 ##########################################################################################
##########################################################################################################

## Definition der Funktion "ModProv" an deren Ende ein vollständiges dist-Objekt steht.
ModProv <- function(df){
  
  ## Zwischenspeichern der übergebenen Markerdaten als Matrix "mat".
  mat <- as.matrix(df)
  
  ## Erstellen der leeren Distanzmatrix
  e <- matrix(0, nrow(df), nrow(df))
  
  ## Zwischenspeichern der Genotypnamen
  genotypes.names <- row.names(df)
  
  ## Zwischenspeichern der Genotypenanzahl
  genotypes.number <- nrow(df)
  
  ## Erstelle Index aller durchzuführenden Genotypenvergleiche ohne Duplikate
  index <- cbind(col(e)[col(e) < row(e)], row(e)[col(e) < row(e)])
  
  ## Erstelle Sequenzen zur Aufteilung der Markerdaten nach Referenz- und SNP-Allel
  ref <- seq(1, ncol(df), 2)
  snp <- seq(2, ncol(df), 2)
  
  ## Definition der Unterfunktion "modprov". Diese Funktion macht die eigentlich Arbeit beim Berechnen der 
  ## Distanzen. Sie wird jeweils für ein Genotypenpaar p,q aufgerufen.
  modprov <- function(x) {
    
    ## Mittels der oben erstellten Sequenzen "ref" und "snp" und dem innerhalb der Funktion als "x" 
    ## übergebenen "index" werden die Markerdaten der Genotypen p und q in
    ## 4 Vektoren gleicher Länge geschrieben, die direkt miteinander verrechnet werden können.
    p.ref <- mat[x[1], ref]
    p.snp <- mat[x[1], snp]
    q.ref <- mat[x[2], ref]
    q.snp <- mat[x[2], snp]
    
    ## "temp" enthält die Ergebnisse aller Markervergleiche zwischen zwei Genotypen. Die Länge der
    ## Vektoren "p.ref", "p.snp", "q.ref" und "q.snp" entspricht der Anzahl der verwendeten Marker.
    ## Für Marker 1 wird die untenstehende Rechnung mit den jeweils an Position 1 der 4 Vektoren
    ## hinterlegten Daten durchgeführt und das Ergebnis an Position 1 von "temp" hinterlegt. Usw. für alle
    ## weiteren Marker.
    temp <- (p.ref + q.ref) * (p.snp + q.snp) * 2^(p.ref * p.snp * q.ref * q.snp - (p.ref + q.snp) * (q.ref + p.snp))
    
    ## Der Mittelwert von "temp" wird gebildet, er stellt die Distanz zwischen Genotyp p und q dar.
    ## Vorhandene Fehlwerte werden zuvor aus dem Vektor "temp" entfernt und nicht berücksichtigt.
    D <- mean(temp, na.rm = TRUE)
    
    ## "D" von Genotyp p und q wird von der Funktion zurückgegeben
    return(D)
  }
  
  ## Fülle die Matrix e mit den genetischen Distanzen zwischen allen Genotypkombinationen indem die Funktion
  ## "modprov" auf alle Zeilen des Objektes "index" angewendet wird (in diesem stehen alle durchzuführenden
  ## Vergleiche, einer pro Zeile)
  e <- unlist(apply(index, 1, modprov))
  
  
  ## Setze einige Attribute für das Objekt "e" um es also dist-Objekt für weitere Funktionen nutzen zu können.
  attr(e, "Size") <- genotypes.number
  attr(e, "Labels") <-  genotypes.names
  attr(e, "Diag") <- FALSE
  attr(e, "Upper") <- FALSE
  attr(e, "method") <- "Modified Prevosti nach Richter/Siebrecht-Schöll"
  attr(e, "call") <- match.call()
  class(e) <- "dist"
  
  ## Gib das fertige "e" als Ergebnis der Funktion aus.
  return(e)
}
############################################################################################################################
