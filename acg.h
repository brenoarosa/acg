/*
 * Programa de analise nodal modificada no regime do tempo e senoidal permanente.
 * Baseado no programa MNA1.c do Professor Antonio Carlos Moreirão de Queiroz (http://www.coe.ufrj.br/~acmq/)
 *
 * Por Breno Arosa - brenoarosa@poli.ufrj.br
 *         Daniel Barradas - danielbarradas@poli.ufrj.br
 *                 Luiz Carlos Macedo - eng.luizmof@poli.ufrj.br
 *                 */

/*
 * Elementos aceitos e linhas do netlist:
 *
 * Resistor:  R<nome> <no+> <no. <resistencia>
 * VCCS:      G<nome> <io+> <io. <vi+> <vi. <transcondutancia>
 * VCVC:      E<nome> <vo+> <vo. <vi+> <vi. <ganho de tensao>
 * CCCS:      F<nome> <io+> <io. <ii+> <ii. <ganho de corrente>
 * CCVS:      H<nome> <vo+> <vo. <ii+> <ii. <transresistencia>
 * Fonte I:   I<nome> <io+> <io. <corrente>
 * Fonte V:   V<nome> <vo+> <vo. <tensao>
 * Amp. op.:  O<nome> <vo1> <vo2> <vi1> <vi2>
 *
 *
 * */

#define versao "1.0 - 17/06/2011"

#define MAX_LINHA			80
#define MAX_NOME			11
#define MAX_ELEM			35
#define MAX_NOS				35
#define TOLG				0.000000000000000000000000000001
#define DEBBUG_ELEMENTOS_
#define DEBBUG_MATRIZ_
#define DEBBUG_NETLIST_
#define DEBBUG_FOURIER_
#define DEBBUG_RESOLVIDO_

#define MAX_TIPO_FONTE		10
#define VALOR_DC			0.00000000000001


#define LAIGUALLB			2
#define SEMLB				3
#define SEMLA				4
#define LANAOINDUTOR		5
#define LBNAOINDUTOR		6

#define DC					1
#define SIN					2
#define PULSADA				3


/*
 * pulso: _/-\
 *
 * tempo 1: _
 * tempo 2: /
 * tempo 3: -
 * tempo 4: \
 * */

using namespace std;

const int MAX_FONTES = 7;
const int MAX_TERMOS = 200;
const int MAX_SOLUCOES = MAX_FONTES *( 2 * MAX_TERMOS + 1 );

typedef long double ld;
typedef int tabela[MAX_NOS+1];

void transcondutancia( complex<ld> gm ,int n1,int n2,int n3,int n4 );
void condutancia( complex<ld> g, int a, int b );
void corrente( complex<ld> i, int a, int b );
void impedancia( complex<ld> valor, int a, int b );
void debbugMatriz (int n);
