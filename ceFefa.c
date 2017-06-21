/*
Circuitos Eletricos 2
Prof. Antonio Carlos Moreirao Queiroz
Grupo: Fernanda F. B. Cassinelli,  Ian Sechin e Rodrigo Ferreira.

Programa de analise de circuitos lineares contendo transistores bipolares,  encontrando o ponto de
operação e a resposta em frequência para pequenos sinais.
*/

/*
Descricao do circuito:
Resistor:      R<nome> <no1> <no2> <resistencia>
Indutor:       L<nome> <no1> <no2> <indutancia>
Acoplamento:   K<nome> <La> <Lb> <k> (La e Lb já declarados)
Capacitor:     C<nome> <no2> <no2> <capacitancia>
VCVC:          E<nome> <noV+> <noV-> <nov+> <nov-> <Av>
CCCS:          F<nome> <noI+> <noI-> <noi+> <noi-> <Ai>
VCCS:          G<nome> <noI+> <noI-> <nov+> <nov-> <Gm>
CCVS:          H<nome> <noV+> <noV-> <noi+> <noi-> <Rm>
Fonte I:       I<nome> <no+> <no-> <modulo> <fase> <valor continuo> (fase em graus)
Fonte V:       V<nome> <no+> <no-> <modulo> <fase> <valor continuo> (fase em graus)
Amp. Op.:      O<nome> <no saida+> <no saida-> <no entrada+> <no entrada->
BJT:           Q<nome> <noc> <nob> <noe> <tipo> <alfa> <alfar> <Isbe> <VTbe> <Isbc> <VTbc> <VA> <C0be> <C1be> <C0bc> <C1bc>
Comentario:    *<comentario>
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_CHAR_LINHA 80
#define MAX_NOME 11
#define MAX_NOME_ARQUIVO 20
#define MAX_ELEMENTOS 500
#define MAX_NOS 50
#define TOLG 1e-30
#define PI 3.14159265358979
#define UM 0.999999999999999999999999999999999999999999
#define ZERO 0.0000000000000000000000000000000000000001


typedef struct infoIndutor {
  double indutancia;
  int id; // indice do elemento na estrutura de netlist
} infoIndutor;

typedef struct infoCapacitor {
  double capacitancia;
  int id;
} infoCapacitor;


typedef struct acoplar {
  char lA[MAX_NOME], lB[MAX_NOME];
} acoplar;

typedef struct infoBJT {
  int noc, nob, noe;
  char tipo[MAX_NOME];
  double alfa, alfar, isbc, vtbc, isbe, vtbe, va, cZerobe, cUmbe, cZerobc,cUmbc;
  int id;
} infoBJT;

typedef struct infoElementos {
  char nome[MAX_NOME], tipo[MAX_NOME], modo[MAX_NOME];
  double valor, modulo, fase;
  int a, b, c, d, x, y;
  double gm,  rm, cZerobe, cUmbe, cZerobc,cUmbc;
  int invertido;
} infoElementos;

acoplar acoplamento[MAX_ELEMENTOS];
infoElementos netlist[MAX_ELEMENTOS];
infoIndutor indutor[MAX_ELEMENTOS];
infoCapacitor capacitor[MAX_ELEMENTOS];
infoBJT bjt[MAX_ELEMENTOS];

int numeroVariaveis, numeroNos, tem, nC, nL, nBJTs, nElementos, listaNos[MAX_NOS];

char 
  escala[3],
  noA[MAX_NOME], noB[MAX_NOME], noC[MAX_NOME], noD[MAX_NOME],
  *p,
  elemento,
  linha[MAX_CHAR_LINHA],
  file[MAX_NOME_ARQUIVO];

double pontos, freqInicial, freqFinal, Yn[MAX_NOS][MAX_NOS], g;
FILE *arquivo;

int no(char *nome){
  int i, existe, numeroNo;
  existe = 0;
  numeroNo = atoi(nome);

  if (numeroVariaveis == 0) {
    listaNos[0] = numeroNo;
    numeroVariaveis++;
    return numeroNo;
  }
  
  for (i = 0; i < numeroVariaveis && existe == 0; i++){
    if (numeroNo == listaNos[i])
      existe = 1;
  }

  if (existe == 1)
    return numeroNo;

  else {
    if (i == MAX_NOS - 1){
      printf("O numero maximo de %u nos foi ultrapassado.\n", MAX_NOS);
      exit (1);
    }
    listaNos[numeroVariaveis] = numeroNo;
    numeroVariaveis++;
  }
  return numeroNo;

}

void montaEstampaDC(void) {
  int i, j;
  char tipo;
  
  for (i = 0; i < numeroVariaveis; i++) {
    for (j=0; j < numeroVariaveis; j++)
      Yn[i][j]=0;
  }
  for (i = 0; i < nElementos; i++) {
        tipo = netlist[i].nome[0];
        if (tipo=='R' || tipo=='C' ) {
          g=1/netlist[i].valor;
          Yn[netlist[i].a][netlist[i].a]+=g;
          Yn[netlist[i].b][netlist[i].b]+=g;
          Yn[netlist[i].a][netlist[i].b]-=g;
          Yn[netlist[i].b][netlist[i].a]-=g;
        }
        else if (tipo=='L'){  //estampa do indutor controlado a corrente (P.O.)
          g=netlist[i].valor;
          Yn[netlist[i].a][netlist[i].x]+=1;
          Yn[netlist[i].b][netlist[i].x]-=1;
          Yn[netlist[i].x][netlist[i].a]-=1;
          Yn[netlist[i].x][netlist[i].b]+=1;
          Yn[netlist[i].x][netlist[i].x]+=g;
        }
        else if (tipo=='G') {
          g=netlist[i].valor;
          Yn[netlist[i].a][netlist[i].c]+=g;
          Yn[netlist[i].b][netlist[i].d]+=g;
          Yn[netlist[i].a][netlist[i].d]-=g;
          Yn[netlist[i].b][netlist[i].c]-=g;
          
        }
        else if (tipo=='I') {
          g=netlist[i].valor;
          Yn[netlist[i].a][numeroVariaveis+1]-=g;
          Yn[netlist[i].b][numeroVariaveis+1]+=g;
        }
        else if (tipo=='V') {
          Yn[netlist[i].a][netlist[i].x]+=1;
          Yn[netlist[i].b][netlist[i].x]-=1;
          Yn[netlist[i].x][netlist[i].a]-=1;
          Yn[netlist[i].x][netlist[i].b]+=1;
          Yn[netlist[i].x][numeroVariaveis+1]-=netlist[i].valor;
        }
        else if (tipo=='E') {
          g=netlist[i].valor;
          Yn[netlist[i].a][netlist[i].x]+=1;
          Yn[netlist[i].b][netlist[i].x]-=1;
          Yn[netlist[i].x][netlist[i].a]-=1;
          Yn[netlist[i].x][netlist[i].b]+=1;
          Yn[netlist[i].x][netlist[i].c]+=g;
          Yn[netlist[i].x][netlist[i].d]-=g;
        }
        else if (tipo=='F') {
          g=netlist[i].valor;
          Yn[netlist[i].a][netlist[i].x]+=g;
          Yn[netlist[i].b][netlist[i].x]-=g;
          Yn[netlist[i].c][netlist[i].x]+=1;
          Yn[netlist[i].d][netlist[i].x]-=1;
          Yn[netlist[i].x][netlist[i].c]-=1;
          Yn[netlist[i].x][netlist[i].d]+=1;
        }
        else if (tipo=='H') {
          g=netlist[i].valor;
          Yn[netlist[i].a][netlist[i].y]+=1;
          Yn[netlist[i].b][netlist[i].y]-=1;
          Yn[netlist[i].c][netlist[i].x]+=1;
          Yn[netlist[i].d][netlist[i].x]-=1;
          Yn[netlist[i].y][netlist[i].a]-=1;
          Yn[netlist[i].y][netlist[i].b]+=1;
          Yn[netlist[i].x][netlist[i].c]-=1;
          Yn[netlist[i].x][netlist[i].d]+=1;
          Yn[netlist[i].y][netlist[i].x]+=g;
        }

        else if (tipo=='Q') {
          g=netlist[i].valor;
          Yn[netlist[i].a][netlist[i].y]+=1;
          Yn[netlist[i].b][netlist[i].y]-=1;
          Yn[netlist[i].c][netlist[i].x]+=1;
          Yn[netlist[i].d][netlist[i].x]-=1;
          Yn[netlist[i].y][netlist[i].a]-=1;
          Yn[netlist[i].y][netlist[i].b]+=1;
          Yn[netlist[i].x][netlist[i].c]-=1;
          Yn[netlist[i].x][netlist[i].d]+=1;
          Yn[netlist[i].y][netlist[i].x]+=g;
        }

    }
  }


int main(void)
{
  printf("Programa de Analise de Ponto de Operacao e Resposta em Frenquencia de Circuitos com BJT\n");
  printf("Por: Fernanda Cassinelli,  Ian Sechin,  Rodrigo Ferreira\n");
  nElementos = 0; nL = 0; nC = 0; nBJTs = 0; numeroVariaveis = 0;
  
  printf("Arquivo com o netlist: ");
  scanf("%50s", file);

  arquivo = fopen(file, "r");
  if (arquivo == 0){
    printf("O arquivo %s nao existe.\n", file);
    exit(1);
  }
  printf("\nLendo o netlist...\n\n");

  fgets(linha, MAX_CHAR_LINHA, arquivo);
  printf("Titulo: %s\n", linha);
  strcpy(listaNos,"0");

  while (fgets(linha, MAX_CHAR_LINHA, arquivo)) {
    nElementos++;
    if (nElementos > MAX_ELEMENTOS) {
      printf("O netlist nao pode ter mais que %d elementos.\n", MAX_ELEMENTOS);
      exit(1);
    }

    linha[0] = toupper(linha[0]);
    elemento = linha[0];
    sscanf(linha, "%10s",  netlist[nElementos - 1].nome); // nome do elemento
    p=linha+strlen(netlist[nElementos - 1].nome);

    if (elemento =='R' || elemento =='L' || elemento =='C') { // R, C, L <nome> <no1> <no2> <valor>
      sscanf(p, "%10s%10s%lg",  noA,  noB, &netlist[nElementos -1].valor);
      if (elemento =='L') {
	      indutor[nL].indutancia = netlist[nElementos -1].valor;
        indutor[nL].id = nElementos -1;
        netlist[nElementos -1].valor = 1e-9; // valor para a analise DC (curto)
        printf("%s %s %s %g\n",  netlist[nElementos -1].nome,  noA,  noB,  indutor[nL].indutancia);
        nL++;
      }
      else if (elemento=='C') {
          capacitor[nC].capacitancia = netlist[nElementos -1].valor;
          capacitor[nC].id = nElementos - 1;
          netlist[nElementos -1].valor = 1e9; // valor para a analise DC (aberto)
          printf("%s %s %s %g\n", netlist[nElementos -1].nome, noA, noB, capacitor[nC].capacitancia);
          nC++;
      }
      else 
         printf("%s %s %s %g\n", netlist[nElementos - 1].nome, noA, noB, netlist[nElementos - 1].valor);

      netlist[nElementos - 1].a = no(noA); // retorna o numero referente a um no ja existente ou um novo no
      netlist[nElementos - 1].b = no(noB);

    }
  	else if (elemento=='I' || elemento=='V'){ // I<nome> <no+> <no-> <modulo> <fase> <valor continuo> (fase em graus)
      sscanf(p, "%10s%10s%lg%lg%lg", noA, noB, &netlist[nElementos - 1].modulo, &netlist[nElementos - 1].fase, &netlist[nElementos - 1].valor);
      printf("%s %s %s %g %g %g\n", netlist[nElementos - 1].nome, noA, noB, netlist[nElementos - 1].modulo, netlist[nElementos - 1].fase, netlist[nElementos - 1].valor);
      netlist[nElementos - 1].a = no(noA);
      netlist[nElementos - 1].b = no(noB);
    }
  
  	else if (elemento=='K') { // K<nome> <La> <Lb> <k> (La e Lb já declarados)
      sscanf(p, "%10s%10s%lg", acoplamento[nElementos - 1].lA, acoplamento[nElementos - 1].lB, &netlist[nElementos - 1].valor);
      printf("%s %s %s %g\n", netlist[nElementos - 1].nome, acoplamento[nElementos - 1].lA, acoplamento[nElementos - 1].lB, netlist[nElementos - 1].valor);
    }
  
    else if (elemento=='G' || elemento=='E' || elemento=='F' || elemento=='H') {
      sscanf(p, "%10s%10s%10s%10s%lg", noA, noB, noC, noD, &netlist[nElementos - 1].valor);
      printf("%s %s %s %s %s %g\n", netlist[nElementos - 1].nome, noA, noB, noC, noD, netlist[nElementos - 1].valor);
      netlist[nElementos - 1].a = no(noA);
      netlist[nElementos - 1].b = no(noB);
      netlist[nElementos - 1].c = no(noC);
      netlist[nElementos - 1].d = no(noD);
    }
    else if (elemento=='O') {
      sscanf(p, "%10s%10s%10s%10s", noA, noB, noC, noD);
      printf("%s %s %s %s %s\n", netlist[nElementos - 1].nome, noA, noB, noC, noD);
      netlist[nElementos - 1].a = no(noA);
      netlist[nElementos - 1].b = no(noB);
      netlist[nElementos - 1].c = no(noC);
      netlist[nElementos - 1].d = no(noD);
    }
    
    else if (elemento=='Q') { // Q<nome> <noc> <nob> <noe> <tipo> <alfa> <alfar> <Isbe> <VTbe> <Isbc> <VTbc> <VA> <C0be> <C1be> <C0bc> <C1bc>
      //srand(time(NULL));
      nBJTs++;
      sscanf(p, "%10s%10s%10s%10s%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg", noA, noB, noC, bjt[nBJTs - 1].tipo, &bjt[nBJTs - 1].alfa, &bjt[nBJTs - 1].alfar, &bjt[nBJTs - 1].isbe, &bjt[nBJTs - 1].vtbe, &bjt[nBJTs - 1].isbc, &bjt[nBJTs - 1].vtbc, &bjt[nBJTs - 1].va, &bjt[nBJTs - 1].cZerobe, &bjt[nBJTs - 1].cUmbe, &bjt[nBJTs - 1].cZerobc, &bjt[nBJTs - 1].cUmbc);
      printf("%10s %10s %10s %10s %10s %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", netlist[nElementos - 1].nome, noA, noB, noC, bjt[nBJTs - 1].tipo, bjt[nBJTs - 1].alfa, bjt[nBJTs - 1].alfar, bjt[nBJTs - 1].isbe, bjt[nBJTs - 1].vtbe, bjt[nBJTs - 1].isbc, bjt[nBJTs - 1].vtbc, bjt[nBJTs - 1].va, bjt[nBJTs - 1].cZerobe, bjt[nBJTs - 1].cUmbe, bjt[nBJTs - 1].cZerobc, bjt[nBJTs - 1].cUmbc);
      strcpy(netlist[nElementos - 1].tipo, bjt[nBJTs - 1].tipo);
      netlist[nElementos - 1].a = no(noA);
      netlist[nElementos - 1].b = no(noB);
      netlist[nElementos - 1].c = no(noC);
    
      bjt[nBJTs - 1].noc = netlist[nElementos - 1].a;
      bjt[nBJTs - 1].nob = netlist[nElementos - 1].b; 
      bjt[nBJTs - 1].noe = netlist[nElementos - 1].c;
      bjt[nBJTs - 1].id = nElementos -1;

      netlist[nElementos - 1].cZerobe = 1e9; // valor para a analise dc
      netlist[nElementos - 1].cUmbe = 1e9;
      netlist[nElementos - 1].cZerobc = 1e9;
      netlist[nElementos - 1].cUmbc = 1e9;
    }
    else if (elemento=='.'){
    	sscanf(p, "%10s %lg %lg %lg", escala, &pontos, &freqInicial, &freqFinal);
    	printf("%s %s %g %g %g", netlist[nElementos - 1].nome, escala, pontos, freqInicial, freqFinal);
    	tem=1;
	  }    
    else if (elemento=='*') { /* Comentario comeca com "*" */
      printf("Comentario: %s", linha);
      nElementos--;
    }
    else {
      printf("Elemento desconhecido: %s\n", linha);
      exit(1);
    }
  }
  fclose(arquivo);
  int i = 0;
  printf("Nos: ");
  for(; i < numeroVariaveis; i++){
    printf("%u", listaNos[i]);
  }

  numeroNos = numeroVariaveis;
  char stringNumero[MAX_NOME], lista[MAX_NOS+1][MAX_NOME+2]; // PRECISA DO +1????????????????????????????????????????????W
  for (i = 0; i < numeroVariaveis; i++){
    sprintf(stringNumero, "%u", listaNos[i]);
    strcpy(lista[i], stringNumero);
  }

  char tipo;

  for (i = 0; i < nElementos; i++) {
    tipo = netlist[i].nome[0];
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O' || tipo=='L') {
      numeroVariaveis++;
      if (numeroVariaveis > MAX_NOS) {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n", MAX_NOS);
        exit(1);
      }
      strcpy(lista[numeroVariaveis - 1],"j"); /* Tem espaco para mais dois caracteres */
      sprintf(stringNumero, "%u", numeroVariaveis - 1);
      strcat(lista[numeroVariaveis - 1], stringNumero);
      netlist[i].x=numeroVariaveis;
    }
    else if (tipo=='H') {
      numeroVariaveis = numeroVariaveis + 2;
      if (numeroVariaveis > MAX_NOS) {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
        exit(1);
      }     
      strcpy(lista[numeroVariaveis - 2],"jx");
      sprintf(stringNumero, "%u", numeroVariaveis - 2);
      strcat(lista[numeroVariaveis - 2], stringNumero);
      netlist[i].x = numeroVariaveis - 2;
      strcpy(lista[numeroVariaveis - 1],"jy");
      sprintf(stringNumero, "%u", numeroVariaveis - 1);
      strcat(lista[numeroVariaveis - 1],stringNumero);
      netlist[i].y = numeroVariaveis - 1;
    }
  }

  /* Lista tudo */
  printf("\nVariaveis internas: \n");
  for (i=0; i < numeroVariaveis; i++)
  printf("%d -> %s\n",i,lista[i]);
  
   /* Monta o sistema nodal modificado */
  if(nBJTs > 0) {
    printf("O circuito e nao linear. Seu modelo linearizado tem %d nos, %d variaveis, %d elementos lineares e %d elementos nao lineares (que se decompoe em %d elementos linearizados)., com nElementos=%d\n",numeroNos,numeroVariaveis,nElementos-nBJTs,nBJTs,nBJTs*7,nBJTs);
  }
  else {
    printf("O circuito e linear. Tem %d nos, %d variaveis e %d elementos\n",numeroNos,numeroVariaveis,nElementos);
  }
  

  return 0;
}