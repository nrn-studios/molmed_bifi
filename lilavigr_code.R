#herunterladen der Vergleichmatrix
load("dynProg.RData")

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
m <- 1
mm <- -1
id <- -2


#Funktion um eine Kostenmatrix K mit bestimmten Werten zu erstellen
K <-
  function(m, mm, id) {
    #Eingabe von den Sequenzen S und T, Match m, Mismatch mm, Insertion und Deletion id
    #Matrix mit Standard Eintrag 0 und Anzahl von Reihen und Spalten entsprechend der Länge der Sequenzen
    M <-
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
        if (basen[i] == "-" & basen[j] == "-") {
          M[i, j] <- 0
        } else {
          #wenn die Basen gleich sind, wird der Wert für Match eingegeben
          if (identical(basen[i], basen[j])) {
            M[i, j] <- m
          } else {
            #wenn an einer Stelle keine Base ist, wird der Wert für Deletion/Insertion eingegeben
            if (basen[i] == "-" | basen[j] == "-") {
              M[i, j] <- id
            } else {
              #wenn die Basen nicht gleich sind, wird der Wert für Mismatch eingegeben
              if (!identical(basen[i], basen[j])) {
                M[i, j] <- mm
              }
            }
          }
        }
      }
    }
    #die Kostenmatrix K wird ausgegeben
    return(M)
  }


#erstellen der Kostenmatrix K mit +1 für Match, -1 für Mismatch und -2 für Insertion/Deletion
K(m, mm, id)


#Funktion für eine Dynamic Programming Matrix (DPM)
DPM <-
  function(R, T, K) {
    #Eingabe der zwei Sequenzen R und T und der Kostenmatrix K
    #Matrix mit Standard Eintrag 0 und Länge und Breite der Länge der Sequenzen
    S <-
      matrix(0,
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
              S[i, j]    <- max(deletion, insertion, match)
            }
          }
        }
      }
    }
    #Traceback
    #Matrix die, die mit den optimalen Sequenzen gefüllt werden soll
    seq <- matrix(0, nrow = 2, ncol = (length(R) + length(T)))
    
    
    h <- length(R) #Zeile von S
    e <- length(T) #Spalte von S
    
    Längste.Sequenz <- 0
    if (h > e) {
      Längste.Sequenz <- h
    } else{
      Längste.Sequenz <- e
    }
    #l läuft die Anzahl der Stellen der längeren Sequenzen von hinten durch
    for (l in Längste.Sequenz:1) {
      if (h == 1 & e == 1) {
        return(seq)
      }
      
      #wenn diagonal zurückgegangen wird und die Differenz einem Match entspricht
      if (S[(h - 1), (e - 1)] == (S[h, e] - m)) {
        #werden beide verglichenen Basen in die Matrix eingesetzt
        seq[1, e] <- R[h]
        seq[2, e] <- T[e]
        e == e - 1
        h == h - 1
      } else {
        #wenn eine Spalte zurückgegangen wird und die Differenz einer Deletion/Insertion entspricht
        if (S[(h - 1), e] == (S[h, e] - id)) {
          #wird in die erste Zeile der Matrix die entsprechende Base der ersten Sequenz eingetragen
          seq[1, e] <- R[h]
          #und in die zweite Zeile ein Leerzeichen
          seq[2, e] <- "-"
          #bei einer Deletion in der zweiten Sequenz darf e nicht weiterlaufen
          e == e - 1
        } else {
          #wenn eine Zeile zurückgegenagen wird und die Differenz einer Deletion/Insertion entspricht
          if (S[h, (e - 1)] == (S[h, e] - id)) {
            #wird in die erste Zeile der Matrix ein Leerzeichen eingegeben
            seq[1, e] <- "-"
            #und in die zweite Zeile die entsprechende base
            seq[2, e] <- T[e]
            #bei einer Deletion in der ersten Sequenz darf h nicht weiterlaufen
            h == h - 1
          }
        }
      }
    }
    #die Matrix mit der optimalen Anordnungen der beiden Sequenzen wird ausgegeben
    return(seq)
  }





DPM(T_1, T_2, K(m, m, id))
str(df)
