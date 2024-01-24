#herunterladen der Vergleichmatrix
load("DynProg.RData")

#mögliche Basen in einem Vektor dargestellt
basen <- c("-", "A", "C", "G", "T")

#die zu testenden Sequenzen zugeordnet und in Vektoren überführt
t_1 <- "-ACGTC"
T_1 <- unlist(strsplit(t_1, split = ""))
t_2 <- "-AGTC"
T_2 <- unlist(strsplit(t_2, split = ""))
s_1 <- "-AGAGTTGCCAAACCCGCT"
S_1 <- unlist(strsplit(s_1, split = ""))
s_2 <- "-AGGGTTGACATCCGTTTT"
S_2 <- unlist(strsplit(s_2, split = ""))

#Variablen für Match, mismatch und Insertion/Deletion
match <- 1
mismatch <- -1
insertion_deletion <- -2


#Funktion um eine Kostenmatrix K mit bestimmten Werten zu erstellen
Kostenmatrix <-
  function(m, mm, id) {
    #Eingabe von den Sequenzen S und T, Match m, Mismatch mm, Insertion und Deletion id
    #Matrix mit Standard Eintrag 0 und Anzahl von Reihen und Spalten entsprechend der Länge der Sequenzen
    Ergebnismatrix <-
      matrix(
        0,
        nrow = (length(basen)),
        ncol = (length(basen)),
        dimnames = list(basen, basen)
      )
    
    #i läuft alle Basen ab und lässt eine gap-spalte
    for (i in 1:length(basen)) {
      #j läuft alle Base ab und lässt eine gap-zeile
      for (j in 1:length(basen)) {
        #wenn keine Base vorhanden ist soll 0 eingegeben werden
        if (basen[i] == "-" &
            basen[j] == "-") {
          Ergebnismatrix[i, j] <- 0
        } else {
          #wenn die Basen gleich sind, wird der Wert für Match eingegeben
          if (identical(basen[i], basen[j])) {
            Ergebnismatrix[i, j] <- m
          } else {
            #wenn an einer Stelle keine Base ist, wird der Wert für Deletion/Insertion eingegeben
            if (basen[i] == "-" |
                basen[j] == "-") {
              Ergebnismatrix[i, j] <- id
            } else {
              #wenn die Basen nicht gleich sind, wird der Wert für Mismatch eingegeben
              if (!identical(basen[i], basen[j])) {
                Ergebnismatrix[i, j] <- mm
              }
            }
          }
        }
      }
    }
    #die Kostenmatrix K wird ausgegeben
    return(Ergebnismatrix)
  }


#erstellen der Kostenmatrix K mit +1 für Match, -1 für Mismatch und -2 für Insertion/Deletion
KM_1 <- Kostenmatrix(match, mismatch, insertion_deletion)
KM_1


#Funktion für eine Dynamic Programming Matrix (DPM)
DPM <-
  function(R, T, K) {
    #Eingabe der zwei Sequenzen R und T und der Kostenmatrix K
    #Matrix mit Standard Eintrag 0 und Länge und Breite der Länge der Sequenzen
    S <- matrix(0,
                nrow = length(R),
                ncol = length(T),
                dimnames = list(R, T))
    
    #i läuft alle Stellen der Sequenz R durch
    for (i in 1:length(R)) {
      #j läuft alle Stellen der Sequenz T durch
      for (j in 1:length(T)) {
        #bei S(0,0) ist die Matrix 0
        if (i == 1 & j == 1) {
          S[i, j] <- 0
        } else {
          #in der ersten Zeile (S(i,0)) wird der Wert aus K mal i genommen
          if (j == 1) {
            S[i, j] <- (i - 1) * (K[R[i], T[j]])
          } else {
            #in der ersten Spalte (S(0,j)) wird der Wert aus K mal j genommen
            if (i == 1) {
              S[i, j] <- (j - 1) * (K[R[i], T[j]])
            } else {
              #in allen anderen Zeilen und Spalten wird der maximale wert aus K von drei Möglichkeiten gesucht (welche passen am besten zsm)
              deletion  <-
                (S[(i - 1), (j)] + K[R[i], "-"]) #die Base eine Spalte davor und die aktuelle Zeile
              insertion <-
                (S[(i), (j - 1)] + K["-", T[j]]) #die Base eine zeile davor und die aktuelle Spalte
              match     <-
                (S[(i - 1), (j - 1)] + K[R[i], T[j]]) #die Basen eine Zeile und eine Spalte davor
              S[i, j]    <-
                max(deletion, insertion, match)
            }
          }
        }
      }
    }
    return(S)
  }

DPM_1 <- DPM(T_1, T_2, KM_1)
DPM_1

DPM_2 <- DPM(T_2, T_1, KM_1)
DPM_2


#Traceback
Traceback <-
  function(S, seq_R, seq_T, K) {
    #Matrix die, die mit den optimalen Sequenzen gefüllt werden soll
    seq <- matrix(0,
                  nrow = 2,
                  ncol = (length(seq_R) + length(seq_T)))
    
    
    h <- length(seq_R) #Zeile von S
    e <- length(seq_T) #Spalte von S
    position <- 1
    
    #solange h und e größer 1 sind, soll die Schleife laufen
    while (h > 0 && e > 0) {
      cat("P:", position, ", h:", h, ", e:", e, "\n")
      if(h == 1 && e == 1 || position == 20){
        break
      }
      
      scoreCurrent <- S[h, e]
      scoreDiag    <- S[h - 1, e - 1]
      scoreUp      <- S[h - 1, e]
      scoreLeft    <- S[h, e - 1]
      
      #wenn diagonal zurückgegangen wird und die Differenz einem Match entspricht
      if (scoreCurrent == (scoreDiag + K[seq_R[h], seq_T[e]]) && h > 1 && e > 1) {
        #werden beide verglichenen Basen in die Matrix eingesetzt
        seq[1, position] <- seq_R[h]
        seq[2, position] <- seq_T[e]
        h <- h - 1
        e <- e - 1
        
        #wenn eine Spalte zurückgegenagen wird und die Differenz einer Deletion/Insertion entspricht
      } else if (scoreCurrent == scoreLeft + K[seq_R[h], "-"] && e > 1) {
        #wird in die erste Zeile der Matrix ein Leerzeichen eingegeben
        seq[1, position] <- seq_R[h]
        #und in die zweite Zeile die entsprechende base
        seq[2, position] <- "-"
        #bei einer Deletion in der ersten Sequenz darf h nicht weiterlaufen
        h <- h - 1
        
        #wenn eine Zeile zurückgegangen wird und die Differenz einer Deletion/Insertion entspricht
      } else { #if (scoreCurrent == scoreUp + K[seq_T[e], "-"] && h > 1) {
        #wird in die erste Zeile der Matrix die entsprechende Base der ersten Sequenz eingetragen
        seq[1, position] <- "-"
        #und in die zweite Zeile ein Leerzeichen
        seq[2, position] <- seq_T[e]
        #bei einer Deletion in der zweiten Sequenz darf e nicht weiterlaufen
        e <- e - 1
      }
      position <- position + 1
    }
    
    #falls eine der beiden Sequenzen vor der anderen endet, muss diese noch aufgefüllt werden
    #also für die erste:
    while (h > 1) {
      seq[1, position] <- seq_R[h]
      seq[2, position] <- "-"
      h <- h - 1
      position <- position + 1
    }
    
    #und für die zweite:
    while (e > 1) {
      seq[1, position] <- "-"
      seq[2, position] <- seq_T[e]
      e <- e - 1
      position <- position + 1
    }
    
    # Entfernen Sie leere Spalten am Ende der Matrix
    seq <- seq[, 1:(position - 1)]
    
    #die Matrix mit der optimalen Anordnungen der beiden Sequenzen wird ausgegeben
    return(seq)
  }

TB_1 <- Traceback(DPM_1, T_1, T_2, KM_1)
TB_1

TB_2 <- Traceback(DPM_2, T_2, T_1, KM_1)
TB_2

str(df)
