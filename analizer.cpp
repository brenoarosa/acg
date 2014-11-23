/*
Programa de analise nodal modificada no regime do tempo e senoidal permanente.
Baseado no programa MNA1.c do Professor Antonio Carlos Moreirão de Queiroz (http://www.coe.ufrj.br/~acmq/)

Por Breno Arosa - brenoarosa@poli.ufrj.br
	Daniel Barradas - danielbarradas@poli.ufrj.br
	Luiz Carlos Macedo - eng.luizmof@poli.ufrj.br
*/

/*
Elementos aceitos e linhas do netlist:

Resistor:  R<nome> <no+> <no. <resistencia>
VCCS:      G<nome> <io+> <io. <vi+> <vi. <transcondutancia>
VCVC:      E<nome> <vo+> <vo. <vi+> <vi. <ganho de tensao>
CCCS:      F<nome> <io+> <io. <ii+> <ii. <ganho de corrente>
CCVS:      H<nome> <vo+> <vo. <ii+> <ii. <transresistencia>
Fonte I:   I<nome> <io+> <io. <corrente>
Fonte V:   V<nome> <vo+> <vo. <tensao>
Amp. op.:  O<nome> <vo1> <vo2> <vi1> <vi2>


*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <complex>
#include "acg.h"

using namespace std;

complex<ld> UM(1.0,0.0);
complex<ld> ZERO(0.0,0.0);
int nTermos; 
ld x;
ld w ;// frequencia
ld respostas[MAX_ELEM];
ld menorPeriodo=9999999;

typedef struct elemento { /* Elemento do netlist */
  char nome[MAX_NOME];
  complex<ld> valor;
  int a,b,c,d,x,y;  // no acoplamento usarei a e b como posicao de La e Lb no netlist
  int tipoFonte;
  ld continuo, amplitude, freq, angulo;
  ld v1,v2,tempoSubida,tempoDescida,tempoLigada,periodo;
} elemento;


typedef struct matrizStruct {
  complex<ld> Yn[MAX_NOS+1][MAX_NOS+2];
  int tipoFonte;
  ld w;
} matrizStruct;


elemento netlist[MAX_ELEM]; /* Netlist */
elemento cacheNetlist[MAX_ELEM];
matrizStruct matriz[MAX_SOLUCOES];

double pi = M_PI;
double PI = M_PI;

int
  n, // Matrizes Yn
  ne, /* Elementos */
  nv, /* Variaveis */
  neq,/* Equacoes */
  nn, /* Nos */
  i,j,k;

tabela C,L;
// C > COLUNAS ; L > LINHAS

char
/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
  nomearquivo[MAX_LINHA+1],
  tipo,
  na[MAX_NOME],nb[MAX_NOME],nc[MAX_NOME],nd[MAX_NOME],
  lista[MAX_NOS+1][MAX_NOME+2], /*Tem que caber jx antes do nome */
  txt[MAX_LINHA+1],
  *p;
FILE *arquivo;

complex<ld> g;

class CFourier{

	ld periodo, tempoLigada, tempoDesligada, tempoSubida, tempoDescida;
	ld Wo;
	ld v1,v2;
	
	ld resultado[MAX_TERMOS][2]; //[tempo][termo][A:0 / B:1]
	
	  
	public:
		CFourier( ld _periodo, ld _tempoDesligada, ld _tempoSubida, ld _tempoLigada, ld _v1, ld _v2);
		int calcularTermos(int termos);
		int printarTermos(int termos);

		ld getWn (int n);
		ld getAo (void);
		ld getAn (int n);
		ld getBn (int n);
		ld getTempoDesligada (void);		
	private:
		void setTempoDesligada (ld periodo,ld tempoSubida, ld tempoLigada, ld tempoDescida);
		void setWo(ld periodo);
		ld setAo (void);
		ld setAn (int n);
		ld setBn (int n);
};

CFourier::CFourier ( ld _periodo,ld _tempoDescida, ld _tempoSubida, ld _tempoLigada, ld _v1, ld _v2){
// Construtor
	periodo = _periodo;
	tempoDescida = _tempoDescida;
	tempoSubida = _tempoSubida;
	tempoLigada = _tempoLigada;
	setTempoDesligada (periodo, tempoSubida, tempoLigada, tempoDescida);
	setWo(periodo);
	v1 = _v1;
	v2 = _v2;
}

void CFourier::setWo (ld periodo){
  Wo = (2*PI)/periodo;
}

ld CFourier::getWn (int n){
  ld result;
  result = n*Wo;
  return result;
}

void CFourier::setTempoDesligada (ld periodo, ld tempoSubida, ld tempoLigada, ld tempoDescida){
	tempoDesligada = periodo - tempoSubida - tempoLigada - tempoDescida;
}

ld CFourier::getTempoDesligada(void){
return tempoDesligada;
}

ld CFourier::getAo(void){
	return resultado[0][0];
}

ld CFourier::getAn(int n){
	return resultado[n][0];
}

ld CFourier::getBn(int n){
	return resultado[n][1];
}

ld CFourier::setAo(void){
  ld result;
  
  result = (v2*tempoLigada - v1*(tempoDescida - periodo + tempoSubida + tempoLigada) + (tempoDescida*(v1 + v2))/2 + (tempoSubida*(v1 + v2))/2)/periodo;
  
  return result;
}

ld CFourier::setAn(int n){
  ld result;
  ld a,b;

  if(tempoSubida==0 && tempoDescida!=0)
	result = (v1*periodo*cos((2*pi*n*(tempoDescida + tempoLigada))/periodo) - v2*periodo*cos((2*pi*n*(tempoDescida + tempoLigada))/periodo) - v1*periodo*cos((2*pi*n*tempoLigada)/periodo) + v2*periodo*cos((2*pi*n*tempoLigada)/periodo) + 2*v1*n*pi*tempoDescida*sin(2*pi*n))/(2*pow(n,2)*pow(pi,2)*tempoDescida);

  else if(tempoSubida!=0 && tempoDescida==0) 
	result = (v1*periodo - v2*periodo - v1*periodo*cos((2*pi*n*tempoSubida)/periodo) + v2*periodo*cos((2*pi*n*tempoSubida)/periodo) - 2*v1*n*pi*tempoSubida*sin((2*pi*n*(tempoSubida + tempoLigada))/periodo) + 2*v2*n*pi*tempoSubida*sin((2*pi*n*(tempoSubida + tempoLigada))/periodo) + 2*v1*n*pi*tempoSubida*sin(2*pi*n))/(2*pow(n,2)*pow(pi,2)*tempoSubida);
	
  else if(tempoSubida==0 && tempoDescida==0)
	result = (v1*sin(2*pi*n) - v1*sin((2*pi*n*tempoLigada)/periodo) + v2*sin((2*pi*n*tempoLigada)/periodo))/(n*pi);

  else
	result = (v1*periodo*tempoDescida - v2*periodo*tempoDescida + v1*periodo*tempoSubida*cos((2*pi*n*(tempoDescida + tempoSubida + tempoLigada))/periodo) - v2*periodo*tempoSubida*cos((2*pi*n*(tempoDescida + tempoSubida + tempoLigada))/periodo) - v1*periodo*tempoSubida*cos((2*pi*n*(tempoSubida + tempoLigada))/periodo) + v2*periodo*tempoSubida*cos((2*pi*n*(tempoSubida + tempoLigada))/periodo) - v1*periodo*tempoDescida*cos((2*pi*n*tempoSubida)/periodo) + v2*periodo*tempoDescida*cos((2*pi*n*tempoSubida)/periodo) + 2*v1*n*pi*tempoDescida*tempoSubida*sin(2*pi*n))/(2*pow(n,2)*pow(pi,2)*tempoDescida*tempoSubida);

  return result;
}

ld CFourier::setBn(int n){
  ld result;
  ld a,b;

  if(tempoSubida==0 && tempoDescida!=0)
	result = (v1*periodo*sin((2*pi*n*(tempoDescida + tempoLigada))/periodo) - v2*periodo*sin((2*pi*n*(tempoDescida + tempoLigada))/periodo) - v1*periodo*sin((2*pi*n*tempoLigada)/periodo) + v2*periodo*sin((2*pi*n*tempoLigada)/periodo) + 2*pi*v2*n*tempoDescida - 2*v1*n*pi*tempoDescida*cos(2*pi*n))/(2*pow(n,2)*pow(pi,2)*tempoDescida);
  
  else if(tempoSubida!=0 && tempoDescida==0)
	result = -(v1*periodo*sin((2*pi*n*tempoSubida)/periodo) - v2*periodo*sin((2*pi*n*tempoSubida)/periodo) - 2*pi*v1*n*tempoSubida - 2*v1*n*pi*tempoSubida*cos((2*pi*n*(tempoSubida + tempoLigada))/periodo) + 2*v2*n*pi*tempoSubida*cos((2*pi*n*(tempoSubida + tempoLigada))/periodo) + 2*v1*n*pi*tempoSubida*cos(2*pi*n))/(2*pow(n,2)*pow(pi,2)*tempoSubida);
	
  else if(tempoSubida==0 && tempoDescida==0)
	result = v2 - v1*cos(2*pi*n) + v1*cos((2*pi*n*tempoLigada)/periodo) - v2*cos((2*pi*n*tempoLigada)/periodo)/(n*pi);

  else
	result =  -(v1*periodo*tempoDescida*sin((2*pi*n*tempoSubida)/periodo) - v2*periodo*tempoDescida*sin((2*pi*n*tempoSubida)/periodo) - v1*periodo*tempoSubida*sin((2*pi*n*(tempoDescida + tempoSubida + tempoLigada))/periodo) + v2*periodo*tempoSubida*sin((2*pi*n*(tempoDescida + tempoSubida + tempoLigada))/periodo) + v1*periodo*tempoSubida*sin((2*pi*n*(tempoSubida + tempoLigada))/periodo) - v2*periodo*tempoSubida*sin((2*pi*n*(tempoSubida + tempoLigada))/periodo) - 2*pi*v1*n*tempoDescida*tempoSubida + 2*v1*n*pi*tempoDescida*tempoSubida*cos(2*pi*n))/(2*pow(n,2)*pow(pi,2)*tempoDescida*tempoSubida);
	
  return result;
}

int CFourier::calcularTermos (int termos){
  int nTermo;

  resultado[0][0] = setAo();
  for ( nTermo = 1; nTermo <= termos ; nTermo++)
  {
    resultado[nTermo][0] = setAn(nTermo);
    resultado[nTermo][1] = setBn(nTermo);
  }

  return 0;
}

int CFourier::printarTermos (int termos){
  int nTermo;

  printf("Periodo:%Lg tempoLigada:%Lg tempoDesligada:%Lg tempoSubida:%Lg tempoDescida:%Lg Wo:%Lg v1:%Lg v2:%Lg \n", periodo, tempoLigada, tempoDesligada, tempoSubida, tempoDescida, Wo, v1,v2);
  printf("Ao: %Lg \n", resultado[0][0]);
  for ( nTermo = 1; nTermo <= termos ; nTermo++)
  {
    printf ("A%d: %Lg \t", nTermo ,resultado[nTermo][0]);
    printf ("B%d: %Lg \n", nTermo ,resultado[nTermo][1]);
  }  
  
  return 0;
}

void numeroComplexo (complex<ld> *x, ld parteReal ,ld parteImag)
{
  complex<ld> numero(parteReal,parteImag);
  *x = numero;
}

complex <ld> numeroComplexo (ld parteReal, ld parteImag)
{
  complex<ld> numero(parteReal,parteImag);
  return numero;
}

void zerarOutrasFontes ( int i)
{
  int n;
  
  for (n=1; n<=ne ; n++)
  {
    if (n != i )
	{
	  tipo = cacheNetlist[n].nome[0];
	  
	  if (tipo == 'V' || tipo == 'I')
	  {
	    cacheNetlist[n].valor = ZERO;
	  }
	}
  }
}

complex <ld> fasor (ld amplitude, ld angulo) // FASOR DE SENO - considerando angulo em radiano
{
  complex <ld> resultado;
  while (angulo > 2*PI) 
	angulo=angulo-2*PI;
  while (angulo < 0 )
	angulo=angulo+2*PI;
	
  resultado = numeroComplexo(amplitude*cos(angulo),amplitude*sin(angulo));
  return resultado;
}

void montarEstampas(elemento netlist[], int n, ld w){
  // Zerar Netlist  
  for (i=0; i<=neq; i++) {
    for (j=0; j<=neq+1; j++)
      matriz[n].Yn[i][j]=0;
  }
    
  /*FAZER UMA FUNÇÃO */
  /* Monta estampas */
   for (i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
    if (tipo=='R'){ 
	condutancia(UM/netlist[i].valor,netlist[i].a,netlist[i].b);
	}
	else if (tipo=='C'){
	  netlist[i].valor *= numeroComplexo(0.0, w);
	  condutancia(netlist[i].valor,netlist[i].a,netlist[i].b);
	}
    else if (tipo == 'L')
    {
	  netlist[i].valor *= numeroComplexo(0.0, w);
	  transcondutancia(UM, netlist[i].a, netlist[i].b, netlist[i].x , 0);
      transcondutancia(UM, 0, netlist[i].x, netlist[i].a, netlist[i].b);
      impedancia(netlist[i].valor,netlist[i].x,netlist[i].x);
	}
	else if (tipo == 'K')
    {
	  netlist[i].valor *= numeroComplexo(0.0, w);
      impedancia(netlist[i].valor, netlist[netlist[i].c].x , netlist[netlist[i].d].x); // lembrando que netlist[ne].c representa a posicao de La no netlist
	  impedancia(netlist[i].valor, netlist[netlist[i].d].x , netlist[netlist[i].c].x);
    }
	
    else if (tipo=='G'){
 	  transcondutancia(netlist[i].valor,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d);
    }
	else if (tipo=='I') corrente(netlist[i].valor,netlist[i].a,netlist[i].b);
    
	else if (tipo=='V') {
      transcondutancia(UM, netlist[i].a, netlist[i].b, netlist[i].x , 0);
      transcondutancia(UM, 0, netlist[i].x, netlist[i].a, netlist[i].b);
      corrente(netlist[i].valor, netlist[i].x ,0);
    }
    else if (tipo=='E') {
      transcondutancia(UM,0,netlist[i].x,netlist[i].a,netlist[i].b);
      transcondutancia(netlist[i].valor,netlist[i].x,0,netlist[i].c,netlist[i].d);
	  transcondutancia(UM, netlist[i].a, netlist[i].b, netlist[i].x, 0);
    }
    else if (tipo=='F') { 
      transcondutancia(netlist[i].valor,netlist[i].a,netlist[i].b,netlist[i].x,0);
      transcondutancia(UM,netlist[i].c,netlist[i].d,netlist[i].x,0);
	  transcondutancia(UM, 0, netlist[i].x, netlist[i].c,netlist[i].d);
    }
    else if (tipo=='H') {
      transcondutancia(UM,0,netlist[i].y,netlist[i].a,netlist[i].b); // L:y ; C:a,b
      transcondutancia(netlist[i].valor,netlist[i].y,0,netlist[i].x,0); // L:y C:x
      transcondutancia(UM,netlist[i].c,netlist[i].d,netlist[i].x,0); // L: c,d C:x
      transcondutancia(UM,0,netlist[i].x,netlist[i].c,netlist[i].d); // L:x ; C:a,b
	  transcondutancia(UM,netlist[i].a,netlist[i].b,netlist[i].y,0); // L:a,b C:y
    }
    else if (tipo=='O') {
	  transcondutancia(UM, netlist[i].x, 0, netlist[i].c, netlist[i].d);
      transcondutancia(UM, netlist[i].a, netlist[i].b, netlist[i].x, 0);
    }
#ifdef DEBUG_ELEMENTOS
	// Printa as Matrizes Finais
    /* Opcional: Mostra o sistema apos a montagem da estampa */
    
	
	printf("Sistema apos a estampa de %s\n",netlist[i].nome);
    debbugMatriz(n);
    getchar();
#endif
  }
}

void copiarNetlist(void){

int j;

// copia netlist para cacheNetlist
for (j=1; j<=ne ; j++){
	strcpy(cacheNetlist[j].nome, netlist[j].nome);
	cacheNetlist[j].valor = netlist[j].valor;  
	cacheNetlist[j].a = netlist[j].a;
	cacheNetlist[j].b = netlist[j].b;
	cacheNetlist[j].c = netlist[j].c;
	cacheNetlist[j].d = netlist[j].d;
	cacheNetlist[j].x = netlist[j].x;
	cacheNetlist[j].y = netlist[j].y;
	cacheNetlist[j].tipoFonte = netlist[j].tipoFonte;
	cacheNetlist[j].continuo = netlist[j].continuo;
	cacheNetlist[j].amplitude = netlist[j].amplitude;
	cacheNetlist[j].freq = netlist[j].freq;
    cacheNetlist[j].angulo = netlist[j].angulo;
	cacheNetlist[j].v1 = netlist[j].v1;
	cacheNetlist[j].v2 = netlist[j].v2;
	cacheNetlist[j].tempoSubida = netlist[j].tempoSubida;
	cacheNetlist[j].tempoDescida = netlist[j].tempoDescida;
	cacheNetlist[j].tempoLigada = netlist[j].tempoLigada;
	cacheNetlist[j].periodo = netlist[j].periodo;	
	}
}

void debbugNetlist(){
  int i;
  
  printf("Netlist parcial da fonte %s\n", cacheNetlist[i].nome);
  for (i=1; i<=ne; i++) {
	tipo=cacheNetlist[i].nome[0];
	if (tipo=='R' || tipo=='C' ) {
	  printf("%s %d %d %+3.4Lg\n",cacheNetlist[i].nome,cacheNetlist[i].a,cacheNetlist[i].b,real(cacheNetlist[i].valor));
	}
	if (tipo=='I' || tipo=='V' ) {
	  printf("%s %d %d %+3.4Lg %+3.4Lgi \n",cacheNetlist[i].nome,cacheNetlist[i].a,cacheNetlist[i].b,real(cacheNetlist[i].valor), imag(cacheNetlist[i].valor));
	}
	else if (tipo == 'L'){
	  printf("%s %d %d %+3.4Lg\n",cacheNetlist[i].nome,cacheNetlist[i].a,cacheNetlist[i].b,real(cacheNetlist[i].valor));
	}
	else if (tipo == 'K'){
	  printf("%s %d %d %+3.4Lg ACOPLAMENTO\n",cacheNetlist[i].nome,cacheNetlist[i].c,cacheNetlist[i].d,real(cacheNetlist[i].valor)); 
	  // AQUI cacheNetlist[I].C E cacheNetlist[I].D REPRESENTAM POSICAO DOS INDUTORES NO cacheNetlist
	}
	else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H') {
	  printf("%s %d %d %d %d %+3.4Lg\n",cacheNetlist[i].nome,cacheNetlist[i].a,cacheNetlist[i].b,cacheNetlist[i].c,cacheNetlist[i].d,real(cacheNetlist[i].valor));
	}
	else if (tipo=='O') {
	  printf("%s %d %d %d %d\n",cacheNetlist[i].nome,cacheNetlist[i].a,cacheNetlist[i].b,cacheNetlist[i].c,cacheNetlist[i].d);
	}
	if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O')
	  printf("Corrente jx: %d\n",cacheNetlist[i].x);
	
	else if (tipo == 'L')
	  printf("Corrente jx: %d\n",cacheNetlist[i].x);
	
	else if (tipo=='H')
	  printf("Correntes jx e jy: %d, %d\n",cacheNetlist[i].x,cacheNetlist[i].y);
	
  }
	  
  getchar();
}

void montarMatrizes(void){
// A partir daqui temos que analisar a série de Fourier  
// A principio a idéia eh pegar cada fonte decompor nos termos , truncar e criar uma matriz de solução para cada um

/* netlist original sera mantido.
   copia-se o netlist para uma variavel local netlist'
   netlist' sera modificada zerando as fontes necessarias
   netlist' sera modificada zerando as fontes necessarias
   monta-se as estampas a partir do netlist'
   matriz[n].Yn[L,C] - n: numero da matriz de estampas atual
   
*/

	int index;
	n = 0;

	for (index = 1;index <= ne; index++){
	  tipo = netlist[index].nome[0];
	  if (((tipo=='V') || (tipo == 'I') ) && (netlist[index].tipoFonte == DC)){
		copiarNetlist();
		w = VALOR_DC;
        cacheNetlist[index].valor = numeroComplexo(cacheNetlist[index].continuo,0.0);
        zerarOutrasFontes(index);					//Zerando outras Fontes
		#ifdef DEBBUG_NETLIST
	        debbugNetlist();
		#endif

		montarEstampas(cacheNetlist, n, w);								//Montando estampa com o resutlado

		#ifdef DEBBUG_MATRIZ
		printf("Matriz n: %d\n", n);
		debbugMatriz(n);
		#endif

		matriz[n].tipoFonte = DC;
		matriz[n].w = w;

		n++;											//Aumentando o índice do matriz.Yn[][]
	  }

	  else if (((tipo=='V') || (tipo == 'I') ) && (netlist[index].tipoFonte == SIN)){
	    copiarNetlist();
		w = VALOR_DC;
		cacheNetlist[index].valor = numeroComplexo(cacheNetlist[index].continuo,0.0);
		zerarOutrasFontes(index);					//Zerando outras Fontes
		#ifdef DEBBUG_NETLIST
			debbugNetlist();
		#endif
		
		montarEstampas(cacheNetlist, n, w);								//Montando estampa com o resutlado
		
		#ifdef DEBBUG_MATRIZ
		printf("Matriz n: %d\n", n);
		debbugMatriz(n);
		#endif
		
		matriz[n].tipoFonte = DC;
		matriz[n].w = w;
		
		n++;
	  
	    copiarNetlist();
        w = cacheNetlist[index].freq * 2 * PI;
		cacheNetlist[index].valor = fasor(cacheNetlist[index].amplitude, cacheNetlist[index].angulo);
        zerarOutrasFontes(index);					//Zerando outras Fontes
        #ifdef DEBBUG_NETLIST
		debbugNetlist();
		#endif
		
		montarEstampas(cacheNetlist, n, w);								//Montando estampa com o resutlado
		
		#ifdef DEBBUG_MATRIZ
		printf("Matriz n: %d\n", n);
		debbugMatriz(n);
		#endif
		
		matriz[n].tipoFonte = SIN;
		matriz[n].w = w;
		
		n++;											//Aumentando o índice do matriz.Yn[][]
	  }
	  
	  else if (((tipo=='V') || (tipo == 'I') ) && (netlist[index].tipoFonte == PULSADA)){
		
		CFourier fou(netlist[index].periodo,
		netlist[index].tempoDescida,
		netlist[index].tempoSubida,
		netlist[index].tempoLigada,
		netlist[index].v1,
		netlist[index].v2);
		
		fou.calcularTermos(nTermos);
		
		// Caso do Ao
		copiarNetlist();
		w = VALOR_DC;
		cacheNetlist[index].valor = numeroComplexo(fou.getAo(),0.0);
		
		#ifdef DEBBUG_FOURIER
		printf(">>>>>>>>>>Termo Ao: %Lg \n", fou.getAo());
		printf(">>>>>>>>>>Fasor: %Lg %Lgi \n", real(cacheNetlist[index].valor), imag(cacheNetlist[index].valor));
		#endif
		
		zerarOutrasFontes(index);					//Zerando outras Fontes
        #ifdef DEBBUG_NETLIST
		debbugNetlist();
		#endif
		
		montarEstampas(cacheNetlist, n, w);								//Montando estampa com o resutlado

		#ifdef DEBBUG_MATRIZ
		printf("Matriz n: %d\n", n);
		debbugMatriz(n);
		#endif
		
		matriz[n].tipoFonte = DC;
		matriz[n].w = w;
		
		n++;	
		
		for (int indexTermos = 1;indexTermos <= nTermos ; indexTermos++){ // loop para varrer todos os termos de An e Bn
			// Caso do An
			copiarNetlist();
			w = fou.getWn(indexTermos);
			cacheNetlist[index].valor = fasor(fou.getAn(indexTermos),PI/2); 
			
			#ifdef DEBBUG_FOURIER
			printf(">>>>>>>>>>Termo A%d: %Lg \n", indexTermos, fou.getAn(indexTermos));
			printf(">>>>>>>>>>Fasor: %Lg %Lgi \n", real(cacheNetlist[index].valor), imag(cacheNetlist[index].valor));
			#endif
			
			zerarOutrasFontes(index);					//Zerando outras Fontes
			#ifdef DEBBUG_NETLIST
			debbugNetlist();
			#endif
			montarEstampas(cacheNetlist, n, w);								//Montando estampa com o resutlado
			
			#ifdef DEBBUG_MATRIZ
			printf("Matriz n: %d\n", n);
			debbugMatriz(n);
			#endif
	
			matriz[n].tipoFonte = SIN;
			matriz[n].w = w;
			
			n++;	
			
			// implementar o Bn ,  similar ao de cima
			copiarNetlist();
			w = fou.getWn(indexTermos);
			cacheNetlist[index].valor = fasor(fou.getBn(indexTermos),0); 
			
			#ifdef DEBBUG_FOURIER
			printf(">>>>>>>>>>Termo B%d: %Lg \n", indexTermos, fou.getBn(indexTermos));
			printf(">>>>>>>>>>Fasor: %Lg %Lgi \n", real(cacheNetlist[index].valor), imag(cacheNetlist[index].valor));
			#endif
			
			zerarOutrasFontes(index);					//Zerando outras Fontes
			#ifdef DEBBUG_NETLIST
			debbugNetlist();
			#endif
			montarEstampas(cacheNetlist, n, w);								//Montando estampa com o resutlado
			
			#ifdef DEBBUG_MATRIZ
			printf("Matriz n: %d\n", n);
			debbugMatriz(n);
			#endif
			
			matriz[n].tipoFonte = SIN;
			matriz[n].w = w;
			
			n++;	
		}
	  
	  
	  }
	
	}
	
	if (n != 0 ) // caso o netlist tenha fonte. no fim de montar estampas tempos que diminui o n de 1 para ficar o numero certo de matrizes
	  n--;
}

int procurarIndutor (char* na, char* nb, int* posNa, int* posNb) 
// retorna indutor ou acoplamento; posNa > posicao de na em netlist; posNb > posicao de nb em netlist
{
  int achouA, achouB;
  
  if (!strcmp(na,nb)) // na = nb
    return LAIGUALLB;
	
  if ((na[0] == 'L') && (nb[0] == 'L'))
  {
    achouA = achouB = 0;
	
	for (i=1; i < ne ;i++) // procura até o ne atual
	{
	  if (!strcmp(netlist[i].nome, na))
	  {
	    *posNa = i;
		achouA =1;
      }
      else if (!strcmp(netlist[i].nome, nb))
	  {
  	    *posNb = i;
			achouB =1;
	  }
	}
	
	if (!achouA)
	  return SEMLA;
	
	if (!achouB)
	  return SEMLB;
	  
	return (1);
  }
  
  if (na[0] == 'L' && nb[0] != 'L')
    return LBNAOINDUTOR;
	
  if (nb[0] == 'L' && na[0] != 'L')
    return LANAOINDUTOR;

  return (0); // nao e acoplamento
}

int acoplamentoMsgError(int err)
{
  switch(err)
  {
    case LAIGUALLB:
	  printf("Acoplamento entre indutores iguais! \n");
	  printf("Erro no netlist! \nPrograma sendo fechado.");
      getchar();
	  return (LAIGUALLB);
	case SEMLB:
	  printf("Indutor Lb nao encontrado! \n");
	  printf("Erro no netlist! \nPrograma sendo fechado.");
      getchar();
	  return (SEMLB);
    case SEMLA:
	  printf("Indutor La nao encontrado! \n");
	  printf("Erro no netlist! \nPrograma sendo fechado.");
      getchar();
	  return (SEMLA);
	case LANAOINDUTOR:
	  printf("La nao indutor! \n");
	  printf("Erro no netlist! \nPrograma sendo fechado.");
      getchar();
	  return (LANAOINDUTOR);
	case LBNAOINDUTOR:
	  printf("Lb nao indutor! \n");
	  printf("Erro no netlist! \nPrograma sendo fechado.");
      getchar();
	  return (LBNAOINDUTOR);  
  }
  
  return 0;  
}
   
void resolverSistema(int n)
{
  int i,j,l, a;
  complex<ld> t, p;

  for (i=1; i<=neq; i++) {
    t=(0.0,0.0);
    a=i;
    for (l=i; l<=neq; l++) {
      if (abs(matriz[n].Yn[l][i])>abs(t)) {
	a=l;
	t=matriz[n].Yn[l][i];
      }
    }
    if (i!=a) {
      for (l=1; l<=neq+1; l++) {
	    p=matriz[n].Yn[i][l];
	    matriz[n].Yn[i][l]=matriz[n].Yn[a][l];
	    matriz[n].Yn[a][l]=p;
      }
    }
    if (abs(t)<TOLG) {
      printf("Sistema singular\n");
      exit(1);
    }
    for (j=neq+1; j>0; j--) {  /* Ponha j>0 em vez de j>i para melhor visualizacao */
      matriz[n].Yn[i][j] /= t;
      p=matriz[n].Yn[i][j];
      for (l=1; l<=neq; l++) {
	    if (l!=i)
	      matriz[n].Yn[l][j]=matriz[n].Yn[l][j]-matriz[n].Yn[l][i]*p;
      }
    }
  }
}

void debbugMatriz (int n)
{
    for (k=1; k<=neq; k++) {
      for (j=1; j<=neq+1; j++)
        if ((real(matriz[n].Yn[k][j])!=0) || (imag(matriz[n].Yn[k][j])!=0) ) printf(" |%+4.1Lg%+4.1Lgi| ",real(matriz[n].Yn[k][j]),imag(matriz[n].Yn[k][j]));
        else printf(" |%+4.1Lg%+4.1Lgi| ", real(ZERO), imag(ZERO));
      printf("\n");
	}
}

void resolverSistemas(void)
{
  int j;
  
  for (i=0; i <= n; i++) // seria <= ??????
  {
    resolverSistema(i);
	
	#ifdef DEBBUG_RESOLVIDO
	printf("Sistema resolvido n:%d \n", i);
	debbugMatriz(i);
	#endif
  }	
}

/* Rotina que conta os nos e atribui numeros a eles */
int numero(char *nome)
{
  int i,achou;

  i=0; achou=0;
  while (!achou && i<=nv)
  {achou=!strcmp(nome,lista[i]);
    if (!(achou)) i++;}
  if (!achou) {
    if (nv==MAX_NOS) {
      printf("O programa so aceita ate %d nos\n",nv);
      exit(1);
    }
    nv++;
    strcpy(lista[nv],nome);
    return nv; /* novo no */
  }
  else {
    return i; /* no ja conhecido */
  }
}

void testarnos(void) {
 if (nv>MAX_NOS) {
   printf("As variaveis extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
   exit(1);
 }
}

void transcondutancia(complex<ld>gm ,int n1,int n2,int n3,int n4) {
  matriz[n].Yn[L[n1]][C[n3]]+=gm;
  matriz[n].Yn[L[n2]][C[n4]]+=gm;
  matriz[n].Yn[L[n1]][C[n4]]-=gm;
  matriz[n].Yn[L[n2]][C[n3]]-=gm;
}

void condutancia(complex<ld> g, int a, int b){
  transcondutancia(g,a,b,a,b);
}

void corrente(complex<ld> i, int a, int b) {
  matriz[n].Yn[L[a]][neq+1]-=i;
  matriz[n].Yn[L[b]][neq+1]+=i;
}

void impedancia(complex<ld> valor, int a, int b){
  matriz[n].Yn[L[a]][C[b]]+=valor;
}

int main(void)
{


  char auxNome[MAX_LINHA + 1];
  char tipoFonte[MAX_TIPO_FONTE]; // lemos tipoFonte como string mas convertemos essa informacao pra int no netlist
  complex<ld> La;
  complex<ld> Lb;
  FILE *arqGrafico;
  
  ld tempo;
  int numTermos=0;
  int nFontes=0;
  
  n = 0;
  
  ld tempoFinal = 10;
  ld passo = 0.001;
  
  char testeLeitura[MAX_LINHA];
  
  system("clear");
  printf("Programa de analise nodal modificada no regime do tempo e permanente senoidal\n");
  printf("Por Breno Arosa, Daniel Barradas e Luiz Carlos Macedo\n");
  for(i=0; i<=MAX_NOS; i++) {C[i]=i; L[i]=i;} /* Inicializa tabelas */
 denovo:
  /* Leitura do netlist */
  ne=0; nv=0; strcpy(lista[0],"0");
  printf("Nome do arquivo com o netlist (ex: mna.net): ");
  scanf("%50s",nomearquivo);
  arquivo=fopen(nomearquivo,"r");
  if (arquivo==0) {
    printf("Arquivo %s inexistente\n",nomearquivo);
    goto denovo;
  }
  
  printf("Lendo netlist:\n");
  fgets(txt,MAX_LINHA,arquivo);
  printf("Titulo: %s",txt);
  while (fgets(txt,MAX_LINHA,arquivo)) {
    ne++; /* Nao usa o netlist[0] */
    if (ne>MAX_ELEM) {
      printf("O programa so aceita ate %d elementos\n",MAX_ELEM);
	  getchar();
      exit(1);
    }
    txt[0]=toupper(txt[0]);
    tipo=txt[0];
    sscanf(txt,"%10s",netlist[ne].nome);
    p=txt+strlen(netlist[ne].nome); /* Inicio dos parametros */
    /* O que e lido depende do tipo */
	
	netlist[ne].tipoFonte = netlist[ne].amplitude = netlist[ne].freq = netlist[ne].angulo = 0;
	netlist[ne].v1=netlist[ne].v2=netlist[ne].tempoSubida=netlist[ne].tempoDescida=netlist[ne].tempoLigada=netlist[ne].periodo=0;
    
	if (tipo =='I' || tipo== 'V'){
	  nFontes++;
	  if(nFontes > MAX_FONTES){
		printf("O netlist excede o numero maximo de fontes %d! \n", MAX_FONTES);
		getchar();
		exit(5);
	  }
	  sscanf(p,"%10s%10s%10s ",na,nb,tipoFonte);
	  printf("%s %s %s %s ",netlist[ne].nome,na,nb, tipoFonte);  

	  p+= (strlen(na) + sizeof(char)) + (strlen(nb) + sizeof(char)) + (strlen(tipoFonte) + sizeof(char)); 

	  // p comeca apontando para o espaco em branco entre o nome na
	  // le na e soma sizeof(char) para compensar o EOS
	  // agora p aponta novamente para o espaco em branco entre na e nb
	  
	  // variavel que aloca termos que nao usamos.
	  int *lixo;
	  lixo = (int*) malloc(sizeof(int));
	  
	  if (!strcmp(tipoFonte,"DC"))
	  { 
	    sscanf(p,"%Lg", &netlist[ne].continuo);
	    netlist[ne].tipoFonte = DC;
		printf("%+3.4Lg", netlist[ne].continuo);
      }	  
	  
	  
	  else if (!strcmp(tipoFonte,"SIN")){
	    netlist[ne].tipoFonte = SIN;
	    sscanf(p, "%Lg %Lg %Lg %d %d %Lg %d",&netlist[ne].continuo, &netlist[ne].amplitude, &netlist[ne].freq, lixo, lixo, &netlist[ne].angulo, lixo);
	    printf("%+3.4Lg %+3.4Lg %+3.4Lg %+3.4Lg",netlist[ne].continuo,netlist[ne].amplitude , netlist[ne].freq, netlist[ne].angulo);
		netlist[ne].angulo *= (PI /180);
	  }     
      
	  else if (!strcmp(tipoFonte,"PULSE")){
		
		netlist[ne].tipoFonte = PULSADA;
	    sscanf(p, "%Lg %Lg %d %Lg %Lg %Lg %Lg %d",&netlist[ne].v1, &netlist[ne].v2, lixo, &netlist[ne].tempoSubida, &netlist[ne].tempoDescida, 
			&netlist[ne].tempoLigada,&netlist[ne].periodo, lixo);
	    printf("%+3.4Lg %+3.4Lg %+3.4Lg %+3.4Lg %+3.4Lg %+3.4Lg",netlist[ne].v1, netlist[ne].v2, netlist[ne].tempoSubida,
		netlist[ne].tempoDescida, netlist[ne].tempoLigada,netlist[ne].periodo);
		if (netlist[ne].tempoSubida==0) netlist[ne].tempoSubida= VALOR_DC;
		if (netlist[ne].tempoDescida==0) netlist[ne].tempoDescida= VALOR_DC;
		
		if (menorPeriodo>netlist[ne].periodo) menorPeriodo=netlist[ne].periodo;
	  
	  }
	  
	  free(lixo);
	  
	  printf("\n");
	  netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }
      	
    else if (tipo=='C'){
      sscanf(p,"%10s%10s%Lg",na,nb,&x);
      netlist[ne].valor=numeroComplexo(x, 0.0);
      printf("%s %s %s %+3.4Lg\n",netlist[ne].nome,na,nb,real(netlist[ne].valor));
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }
	
	else if (tipo == 'K'){
	  sscanf(p,"%10s%10s%Lg",na,nb,&x);

	  if (procurarIndutor(na, nb, &netlist[ne].c, &netlist[ne].d) == 1){
		La = netlist[netlist[ne].c].valor; // lembrando que netlist[ne].c representa a posicao de La no netlist
		Lb = netlist[netlist[ne].d].valor;
		netlist[ne].valor = numeroComplexo( ( x*sqrt(real(La))*sqrt(real(Lb)) ), 0.0 ) ;
	    printf("%s %s %s k:%Lg M:%+3.4Lg ACOPLAMENTO\n",netlist[ne].nome,na,nb, x, real(netlist[ne].valor));
	  }
	  
	  else {
	    exit( acoplamentoMsgError( procurarIndutor(na, nb, &netlist[ne].c, &netlist[ne].d) ) );
	  }
    }
	
    else if (tipo == 'L'){ 
      sscanf(p,"%10s%10s%Lg",na,nb,&x);
	  
	  netlist[ne].valor=numeroComplexo(x, 0.0);
	  printf("%s %s %s %+3.4Lg\n",netlist[ne].nome,na,nb,real(netlist[ne].valor));
	  netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);	  
    }

    else if (tipo=='R') {
      sscanf(p,"%10s%10s%Lg",na,nb,&x);
      netlist[ne].valor=x; // está passando somente o valor da parte REAL
      printf("%s %s %s %+3.4Lg\n",netlist[ne].nome,na,nb,real(netlist[ne].valor));
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }
    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H') {
      sscanf(p,"%10s%10s%10s%10s%Lg",na,nb,nc,nd,&x);
	  netlist[ne].valor = numeroComplexo( x, 0.0);
      printf("%s %s %s %s %s %+3.4Lg\n",netlist[ne].nome,na,nb,nc,nd,real(netlist[ne].valor));
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      netlist[ne].c=numero(nc);
      netlist[ne].d=numero(nd);
    }
    else if (tipo=='O') {
      sscanf(p,"%10s%10s%10s%10s",na,nb,nc,nd);
      printf("%s %s %s %s %s\n",netlist[ne].nome,na,nb,nc,nd);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      netlist[ne].c=numero(nc);
      netlist[ne].d=numero(nd);
    }
    else if (tipo=='*') { /* Comentario comeca com "*" */
      printf("Comentario: %s",txt);
      ne--;
    }
	else if (!strcmp(netlist[ne].nome, ".TRAN") ){
	  sscanf(p,"%Lg%Lg%i", &tempoFinal, &passo, &numTermos);
	  ne--;
	}
    else {
      printf("Elemento desconhecido: %s\n",txt);
      getchar();
      exit(1);
    }
  }
  fclose(arquivo);

  int limiteTermos;

  limiteTermos = floor( menorPeriodo/(4*PI*passo) );
  
  if (limiteTermos > MAX_TERMOS)
    limiteTermos = MAX_TERMOS;
	
  if ( numTermos == 0 )
    nTermos = limiteTermos;
  
  else if (numTermos > limiteTermos){
    nTermos = limiteTermos;
	printf("O numero de Termos da Serie de F. ultrapassa o que se pode representar com o passo dado.\n");
  }
  
  else
	nTermos=numTermos;

  if (nTermos > MAX_TERMOS){
	nTermos=MAX_TERMOS;
	printf("O numero de Termos da Serie de F. estourou o limite!\n");
  }
  printf("Utilizando numero de Termos da Serie de F. %i \n",nTermos);
  
  
  
  nn=nv; 
  neq=nn;
  
  /* Atualiza as tabelas e acrescenta variaveis de corrente acima dos nos, anotando no netlist */
  for (i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
    if (tipo=='V') {
      nv++;
      testarnos();
      strcpy(lista[nv],"j"); /* Tem espaco para mais dois caracteres */
      strcat(lista[nv],netlist[i].nome);
      netlist[i].x=nv;
    }
    else if (tipo == 'L')
	{
	  nv++;
      testarnos();
      strcpy(lista[nv],"j"); /* Tem espaco para mais dois caracteres */
      strcat(lista[nv],netlist[i].nome);
      netlist[i].x=nv;
    }

    else if (tipo=='O') {
	  nv++;
      testarnos();
	  strcpy(lista[nv],"j"); /* Tem espaco para mais dois caracteres */
      strcat(lista[nv],netlist[i].nome);
	  netlist[i].x=nv;
    }
	
    else if (tipo=='E') {
      nv++;
      testarnos();
      strcpy(lista[nv],"j"); /* Tem espaco para mais dois caracteres */
      strcat(lista[nv],netlist[i].nome);
      netlist[i].x=nv;
    }
    else if (tipo=='F') {
      nv++;
      testarnos();
      strcpy(lista[nv],"j"); /* Tem espaco para mais dois caracteres */
      strcat(lista[nv],netlist[i].nome);
      netlist[i].x=nv;
    }
	
    else if (tipo=='H') {
      nv=nv+2;
      testarnos();
      strcpy(lista[nv-1],"jx"); strcat(lista[nv-1],netlist[i].nome);
      netlist[i].x=nv-1;
      strcpy(lista[nv],"jy"); strcat(lista[nv],netlist[i].nome);
      netlist[i].y=nv;
    }
  }
  getchar();
  /* Lista tudo */
  printf("Variaveis internas: \n");
  for (i=0; i<=nv; i++)
    printf("%d . %s (%d)\n",i,lista[i],C[i]);
  getchar();
  printf("Netlist interno final\n");
  for (i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
	
	if(tipo=='I' || tipo=='V'){
	  if (netlist[i].tipoFonte == DC) 
	    printf("%s %d %d DC %+3.4Lg\n",netlist[i].nome,netlist[i].a,netlist[i].b, netlist[i].continuo);
	
	  else if (netlist[i].tipoFonte == SIN)
	    printf("%s %d %d SIN %+3.4Lg %+3.4Lg %+3.4Lg %+3.4Lg\n",netlist[i].nome, netlist[i].a, netlist[i].b, netlist[i].continuo , netlist[i].amplitude , netlist[i].freq , netlist[i].angulo);
		
	  else if ( netlist[i].tipoFonte == PULSADA)
	    printf("%s %d %d PULSADA %+3.4Lg %+3.4Lg %+3.4Lg %+3.4Lg %+3.4Lg %+3.4Lg\n",netlist[i].nome, netlist[i].a, netlist[i].b,
			netlist[i].v1, netlist[i].v2, netlist[i].tempoSubida,netlist[i].tempoDescida,netlist[i].tempoLigada, netlist[i].periodo);
	}
	
    if (tipo=='R' || tipo=='C' ) {
      printf("%s %d %d %+3.4Lg\n",netlist[i].nome,netlist[i].a,netlist[i].b,real(netlist[i].valor));
    }
    else if (tipo == 'L'){
      printf("%s %d %d %+3.4Lg\n",netlist[i].nome,netlist[i].a,netlist[i].b,real(netlist[i].valor));
    }
	else if (tipo == 'K'){
	  printf("%s %d %d %+3.4Lg ACOPLAMENTO\n",netlist[i].nome,netlist[i].c,netlist[i].d,real(netlist[i].valor)); 
	  // AQUI NETLIST[I].C E NETLIST[I].D REPRESENTAM POSICAO DOS INDUTORES NO NETLIST
	}
    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H') {
      printf("%s %d %d %d %d %+3.4Lg\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d,real(netlist[i].valor));
    }
    else if (tipo=='O') {
      printf("%s %d %d %d %d\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d);
    }
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O')
      printf("Corrente jx: %d\n",netlist[i].x);
    
    else if (tipo == 'L')
      printf("Corrente jx: %d\n",netlist[i].x);
    
    else if (tipo=='H')
      printf("Correntes jx e jy: %d, %d\n",netlist[i].x,netlist[i].y);
	
  }
  
  if (!strcmp(netlist[ne+1].nome, ".TRAN"))
	  printf("Tempo de simulacao: 0 a %Lg (s) \t passo: %Lg (s) \n", tempoFinal, passo);
	  
  getchar();
  
  neq = nv;
  
  printf("O circuito tem %d nos, %d variaveis, %d equacoes e %d elementos\n",nn,nv,neq,ne);
  getchar();

  montarMatrizes();
  resolverSistemas();

  // Comecando a escrever arquivo!
  strncpy(auxNome, nomearquivo,(strlen(nomearquivo)-4));
  auxNome[strlen(nomearquivo)-4] = '\0';
  strcat(auxNome, ".tab");
	
  printf("Criando arquivo com as solucoes: %s \n",auxNome);
  
  arqGrafico=fopen(auxNome,"w");
  
  // printando primeira linha
  fprintf(arqGrafico, "t ");
  for (i = 1; i <= nv ; i++){
    fprintf(arqGrafico, "%s ", lista[i]);
  }
  fprintf(arqGrafico, "\n");
  

  for (tempo= 0.0; tempo < tempoFinal ; tempo+=passo){	//fazendo o tempo variar , e para cada tempo tenho que calcular as respostas de cada variavel
	fprintf(arqGrafico,"%+3.4Lg ",tempo);	

	for (int var=0; var<=nv; var++)//zerando o array em que eu salvo as respostas
	  respostas[var]=0;
	  
	  
	for (int solucaoN = 0 ; solucaoN <= n ; solucaoN++){
	
		if (matriz[solucaoN].tipoFonte == DC){ // inclui parte DC de fontes senoidais
		  for (int var=1; var <= neq ; var++) {
		    respostas[var] +=  real(matriz[solucaoN].Yn[var][neq+1]);
		  }
		}
		
		if (matriz[solucaoN].tipoFonte == SIN)
		{
			for (int var=1; var <= nv; var++) {
				respostas[var] += ( abs(matriz[solucaoN].Yn[var][neq+1]) * sin( matriz[solucaoN].w * tempo + arg(matriz[solucaoN].Yn[var][neq+1]) ) );
				// resposta += A*sin(wt + teta)
			}
			
		}
		
	}
	
	for (int var=1; var <= nv; var++) {
		fprintf(arqGrafico,"%+9.9Lg ",respostas[var]);
	}
	
	fprintf(arqGrafico,"\n");
  }

  fclose(arqGrafico);
  printf("Arquivo criado! \n");
  getchar();
  return 0;
}
