\\ **********************************************************************************
\\  Fonctions de l'exercice :                                                       *
\\                                                                                  *
\\    encodegln(s,n): Encodage dans Gl(n) de la chaine de caractère s               *
\\    decodegln(M): Décodage de la matrice M dans Gl(n)                             *
\\    ordre(M, p): ordre de la matrice M dans le groupe mult. (Z/pZ)                *
\\ **********************************************************************************

encodegln(s,n)={
  my(v);
  v=[if(x==32,0,x-96)|x<-Vec(Vecsmall(s))];
  if(#v>n^2,warning("string truncated to length ",n^2));
  v = Vec(v,n^2);
  matrix(n,n,i,j,v[(i-1)*n+j]);
}

decodegln(M)={
  my(v);
  v= concat(Vec(M~))~;
  v= Strchr([ if (s == 0, 32, s + 96) | s <- v]);
  v;
}


ordre(M, p) = {
    my(k, n, tmp, M_id_mod);
    k= 1;
    n= matsize(M)[1];

    tmp= Mod(M, p);
    M_id_mod= Mod(matid(n), p);
    while( M_id_mod != tmp, tmp*= M; k++);
    k;
}

\\ **********************************************************************************
\\  Paramètres de l'exercice :                                                      *
\\                                                                                  *
\\    C: matrice d'ordre k*k chiffré du message <clair> à déchiffrer                *
\\    n: exposant de la puissance effectué utilisé pour chiffrer le message <clair> *
\\    k: ordre de la matrice C                                                      *
\\    p: Espace fini dans lequel on chiffre notre message  (Z/pZ)                   *
\\    d: Ordre de la matrice C dans Z/pZ                                            *
\\    u,v: Coefficient de Bezout retrouvé par l'algorithme d'Euclide étendu         *
\\    m: message <clair> au format string                                           *
\\ **********************************************************************************


\\ Initialisation des paramètres requis

text=readstr("input.txt")[1];
n= 65537; \\ 2^16 + 1
k= round(sqrt(#text));
p= 27;
C= encodegln(text, k);


\\ Étape 1 - Trouver l'ordre d de C dans Z/pZ afin d'obtenir M^d = 1 mod p
d= ordre(C, p);

\\ Désormais, nous savons que M^d = 1 mod p

\\ Étape 2 - Effectuer l'algorithme d'Euclide étendu sur l'exposant n et l'ordre d
[u, v]= gcdext(n, d);

\\ **** Nous avons u*n + v*d = 1
\\ **** En l'occurence : u*n = 1 - v*d = 1
\\ **** Donc u est l'inverse de n modulo l'ordre d.
\\ **** On en déduit que : 
\\ ****** C^u = M^(u*n) = M * ( M^(vd) ) = M * ( M^(d) ) ^ v = M * id(k) = M mod p
\\ ****** Avec id(k) la matrice identité d'ordre k

\\ **** C^u = M mod p, il nous suffit d'élever à la puissance u notre matrice de chiffré dans (Z/pZ)

\\ Étape 3 - Elevation à la puissance et conversion (Z/pZ) <---> str


C= lift( Mod(C, 27)^u ); 
m= decodegln(C);

print(m);