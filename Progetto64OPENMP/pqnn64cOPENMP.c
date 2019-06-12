/**************************************************************************************
 * 
 * CdL Magistrale in Ingegneria Informatica
 * Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2018/19
 * 
 * Progetto dell'algoritmo di Product Quantization for Nearest Neighbor Search
 * in linguaggio assembly x86-32 + SSE
 * 
 * Fabrizio Angiulli, aprile 2019
 * 
 **************************************************************************************/

/*
 
 Software necessario per l'esecuzione:

     NASM (www.nasm.us)
     GCC (gcc.gnu.org)

 entrambi sono disponibili come pacchetti software 
 installabili mediante il packaging tool del sistema 
 operativo; per esempio, su Ubuntu, mediante i comandi:

     sudo apt-get install nasm
     sudo apt-get install gcc

 potrebbe essere necessario installare le seguenti librerie:

     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
     sudo apt-get install libc6-dev-i386

 Per generare il file eseguibile:

 nasm -f ef32 pqnn32.nasm && gcc -O0 -m32 -msse pqnn32.o pqnn32c.c -o pqnn32c && ./pqnn32c
 
 oppure
 
 ./runpqnn32

*/

#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>
#include <stdbool.h>

#define 	MATRIX		float*
#define 	x_query 	input->qs
#define     nodo		(input->m+1)
#define 	BLOCKSIZE 	32
#define		p			8
#define		unroll		4


typedef struct {

	char* filename; //
	MATRIX ds; // data set 
	MATRIX qs; // query set
	int n; // numero di punti del data set
	int d; // numero di dimensioni del data/query set
	int nq; // numero di punti del query set
	int knn; // numero di ANN approssimati da restituire per ogni query
	int m; // numero di gruppi del quantizzatore prodotto
	int k; // numero di centroidi di ogni sotto-quantizzatore
	int kc; // numero di centroidi del quantizzatore coarse
	int w; // numero di centroidi del quantizzatore coarse da selezionare per la ricerca non esaustiva
	int nr; // dimensione del campione dei residui nel caso di ricerca non esaustiva
	float eps; // 
	int tmin; //
	int tmax; //
	int exaustive; // tipo di ricerca: (0=)non esaustiva o (1=)esaustiva
	int symmetric; // tipo di distanza: (0=)asimmetrica ADC o (1=)simmetrica SDC
	int silent;
	int display;
	int sub; //numero di dimensioni dei sotto gruppetti (d/m)
	// nns: matrice row major order di interi a 32 bit utilizzata per memorizzare gli ANN
	// sulla riga i-esima si trovano gli ID (a partire da 0) degli ANN della query i-esima
	//probabilment esono long
	int* ANN; 
} params;

int size=p*unroll;
int c_max_heap=0;
float pre_max_heap=0;
bool  nmod4=false;
bool  dmod4=false;
bool  nmod4noex=false;
bool submod4=false;
float max_f=FLT_MAX;

//variabili utili per l'algoritmo esaustivo
float * centroids;
int * pq; 
float* uj_x;
int* c_x;


//variabili utili per l'algoritmo non esaustivo
MATRIX Cc;
int* Cc_index;
float* Cp;
int * Cp_index;
float* stored_distance;
int* len_IL;
int ** IL;
int* bucket;




/*
 * 
 *	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
 * 	mediante un array (float*), in modo da occupare un unico blocco
 * 	di memoria, ma a scelta del candidato possono essere 
 * 	memorizzate mediante array di array (float**).
 * 
 * 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
 * 	matrici per righe (row-major order) o per colonne (column major-order).
 *
 * 	L'assunzione corrente è che le matrici siano in row-major order.
 * 
 */


void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,32); 
}


void free_block(void* p1) { 
	_mm_free(p1);
}


MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(float),rows*cols);
}

int* alloc_vector(int index){
	return (int*) get_block(sizeof(int),index);
}

void dealloc_vector(int* vector) {
	free_block(vector);
}


void dealloc_matrix(MATRIX mat) {
	free_block(mat);
}

/*
 * 
 * 	load_data
 * 	=========
 * 
 *	Legge da file una matrice di N righe
 * 	e M colonne e la memorizza in un array lineare in row-major order
 * 
 * 	Codifica del file:
 * 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
 * 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
 * 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
 * 
 *****************************************************************************
 *	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
 * 	della matrice. 
 *****************************************************************************
 * 
 */

//salvataggio matrice per colonne personalizzato
MATRIX load_data_col_p(char* filename, int *n, int *d, int nn, int dd) {	
	FILE* fp;
	int rows, cols, status, i, cnt;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL) {
		printf("'%s' : bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	rows = nn;
	cols = dd;
		
	MATRIX data1 = alloc_matrix(rows,cols);
	status = fread(data1, sizeof(float), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*d = cols;

	MATRIX data = alloc_matrix(rows,cols);
	for(int j=0; j<cols; j++){
		for(int z=0; z<rows; z++){
			data[z+j*rows] = data1[z*cols+j];
		}
	}
	return data;
}


//salvataggio matrice per colonne
MATRIX load_data_col(char* filename, int *n, int *d) {	
	FILE* fp;
	int rows, cols, status, i, cnt;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL) {
		printf("'%s' : bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
		
	MATRIX data1 = alloc_matrix(rows,cols);
	status = fread(data1, sizeof(float), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*d = cols;

	MATRIX data = alloc_matrix(rows,cols);
	for(int j=0; j<cols; j++){
		for(int z=0; z<rows; z++){
			data[z+j*rows] = data1[z*cols+j];
		}
	}
	return data;
}

MATRIX load_data_row_p(char* filename, int *n, int *d, int nn, int dd) {	
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL) {
		printf("'%s' : bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	rows = nn;
	cols = dd; 
		
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(float), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*d = cols;

	return data;
}

MATRIX load_data_row(char* filename, int *n, int *d) {	
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL) {
		printf("'%s' : bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
		
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(float), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*d = cols;

	return data;
}



void save_ANN(char* filename, int* ANN, int nq, int knn) {	
	FILE* fp;
	int i, j;
	char fpath[256];
	
	sprintf(fpath, "%s.ann", filename);
	fp = fopen(fpath, "w");
	for (i = 0; i < nq; i++) {
		for (j = 0; j < knn; j++)
			fprintf(fp, " %d  ", ANN[i*knn+j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void printX(MATRIX x,int i,int d){
	int j;
	for(j=0;j<d;j++){
		printf("X[%d][%d] = %f\n",i,j,x[i*d+j]);
	}
}

//metodo che stampa i centroidi e i punti associati ad esso (sotto forma di vettore di indici)
void printCentroids(MATRIX C,int* index, int n,int d,int k){
	int i,j,jj;
	for (i=0;i<k;i++){
		for(j=0;j<d;j++){
			printf("C[%d][%d]=   %f \n",i,j,C[i*d+j]);
		}
		printf("il centroide num: %d è associato ai punti --> [",i);
		for(jj=0; jj < n; jj++)
		{
			if (index[jj]==i) {
				printf(".%d. ",jj);
			}
		}
		printf("]\n");
	}
}
//print per stampare il dataset ed il query set
void printDsQs(MATRIX ds, MATRIX qs,int n,int d,int nq){
	printf("n=%d------d=%d---------nq=%d\n",n,d,nq);
	int i,j;
	printf("-------------Data Set-------------\n");
	for (i=0;i<n;i++){
		for(j=0;j<d;j++){
			//printf("ds[%d][%d]=   %f   \n",i,j,ds[i*d+j]); 					//per riga
			printf("ds[%d][%d]=   %f   \n",i,j,ds[i+j*n]);				//per colonna
		}
	}
	if(qs!=NULL){
		printf("-------------Query Set-------------\n");
		for (i=0;i<nq;i++){
			for(j=0;j<d;j++){
				printf("qs[%d][%d]=   %f   \n",i,j,qs[i*d+j]);
			}
		}
	}
}

void printVector(int * v,int n){
	int i;
	for(i=0;i<n;i++){
		printf("v[%d] = %d\n",i,v[i]);
	}
	printf("\n");
}

void printVectorfloat(float * v,int n){
	int i;
	for(i=0;i<n;i++){
		printf("\nv[%d] = %f",i,v[i]);
	}
	printf("\n");
}

//print per testare il metodo Uj
void printEq_row(MATRIX m1, MATRIX m2, int m1_n,int m1_d,int m2_n,int m2_d){
	printf("Si ipotizza che m1 sia più grande di m2 (sia d che n)");
	int i,j;
	for (i=0;i<m1_n;i++){
		for(j=0;j<m1_d;j++){
			
			if(i < m2_n && j < m2_d )
				printf("m1[%d][%d]=   %f   ------>   m2[%d][%d]=  %f  \n",i,j,m1[i*m1_d+j],i,j,m2[i*m2_d+j]);
			else
			{
				printf("m1[%d][%d]=   %f   ------>   m2[%d][%d]= Null \n",i,j,m1[i*m1_d+j],i,j);
			}
			
		}
	}
}

//print per testare il metodo Uj
void printEq_col(MATRIX m1, MATRIX m2, int m1_n,int m1_d,int m2_n,int m2_d){
	printf("Si ipotizza che m1 sia più grande di m2 (sia d che n)");
	int i,j;
	for (i=0;i<m1_n;i++){
		for(j=0;j<m1_d;j++){
			
			if(i < m2_n && j < m2_d )
				printf("m1[%d][%d]=   %f   ------>   m2[%d][%d]=  %f  \n",i,j,m1[i+m1_n*j],i,j,m2[i+m2_n*j]);
			else
			{
				printf("m1[%d][%d]=   %f   ------>   m2[%d][%d]= Null \n",i,j,m1[i+m1_n*j],i,j);
			}
			
		}
	}
}


/* metodi assembly per il precalcolo delle distanze in sdc e adc splittati in Aligned and Unaligned,utilizzati in pre_sdc/pre_adc etc..*/
extern void rowDistance64AdcA(float* c,float* x,float* distance,int i,int j,int k,int sub);
extern void rowDistance64AdcU(float* c,float* x,float* distance,int i,int j,int k,int sub);
extern void rowDistance64SdcA(float* c,float* distance,int i,int j,int j_d,int k,int sub);
extern void rowDistance64SdcU(float* c,float* distance,int i,int j,int j_d,int k,int sub);
/* metodi utilizzati nel k_means per il calcolo delle distanze con i centroidi, presupponendo che il dataset sia salvato per colonna*/
extern void colDistance64Sing(float * ds,float* c,float* distance,int i,int j,int d,int n);
extern void colDistance64U(float * ds,float* c,float* distance,int i,int j,int d,int n);
extern void colDistance64A(float * ds,float* c,float* distance,int i,int j,int d,int n);
//extern void colDistance64Block(float * ds,float* c,float* distance,int i,int j,int d,int n,int b);
extern void	colDistance64OptimizedU(float* data,float* centroids,float* distance,int i,int j,int d,int n);
extern void	colDistance64OptimizedA(float* data,float* centroids,float* distance,int i,int j,int d,int n);
extern void distanceControl64Sing(float * distance,float * min_distance,int *labels,int j,int i);
extern void distanceControl64(float * distance,float * min_distance,int *labels,int j,int i);
//extern void distanceControl64Block(float * distance,float * min_distance,int *labels,int j,int i);

/*metodi utilizzati in k_means per aggiornare i centroidi,ripulire i centroidi ed assegnare al vettore min_distance il valore FLT_MAX*/
extern void updateCentroid(float* c,float* c1,float* counts,int k,int d);
extern void clearCentroids(float* counts,float* c1,int k,int d);
extern void assignValue(float* list,float* value,int i);

/*metodo utilizzato in extract_col per prelevare un sottoinsieme di valori da dataset, presupponendo il salvataggio per colnna*/
extern void extr_col(float* ds, int n, int d, int nr, int divi, float* result);

/*metodi utilizzati per il calcolo delle distanze tra 2 punti,presupponendo il salvataggio per riga*/
extern void dist64A(float * x,float * y,float* distance, int d);
extern void dist64U(float * x,float * y,float* distance, int d);

/*metodi per calcolare il centroide più vicino ad una query x*/
extern void cent_XA(float* cent, float* xx, int k, int dd, float* tmp, int* park, float* dis);
extern void cent_XU(float* cent, float* xx, int k, int dd, float* tmp, int* park, float* dis);

//extern void mapping64(int i, int j, int n, int * indice, int index);


/*metodo per estrarre dal dataset un sottoinsieme di nr righe, presupponendo il salvataggio per riga*/
MATRIX extrac_row(MATRIX ds,int n,int d,int nr){																	
	MATRIX result = alloc_matrix(nr,d);
	int i,j,h;
	for (h = i = 0; i < nr; h += n/nr , i++) {
		for (j = 0; j < d;j++){
			result[i*d+j] = ds[h*d+j];
		}
	}
	return result;
}
/*metodo per estrarre dal dataset un sottoinsieme di nr righe, presupponendo il salvataggio per colonna*/
MATRIX extrac_col(MATRIX ds,int n,int d,int nr){																	
	MATRIX result = alloc_matrix(nr,d);
	int divi = n/nr;
	int i,j,h;
	for (j = 0; j < d;j++) {

		for (i = 0; i < nr; i++){
		
			result[i+nr*j] = ds[i+n*j];
		}
	}
	return result;
}

/**
 * metodo che prende un sottogruppo (sub dimensionale) del data set
 * j serve per prendere il j-esimo gruppetto di dimensioni j=2 equivale ad U2
 * presupponendo ds salvato per colonna*/
MATRIX Uj_col(MATRIX ds, int j,int m,int n,int d){
	int j_d,k,c,i;
	int sub=d/m;
	MATRIX uj = alloc_matrix(n,sub);
	c=0;																					
	for(j_d=j*sub; j_d<(j+1)*sub; j_d++){																	
		for(i=0; i<n; i++){
			uj[i+n*c] = ds[i+n*j_d];
			
		}
		c++;
	}																					
	return uj;
}
/**
 * metodo che prende un sottogruppo (sub dimensionale) del data set
 * j serve per prendere il j-esimo gruppetto di dimensioni j=2 equivale ad U2
 * presupponendo ds salvato per riga*/
MATRIX Uj_row(MATRIX ds, int j,int m,int n,int d){
	int j_d,k,c,i;
	int sub=d/m;
	MATRIX uj = alloc_matrix(n,sub);
	for(i = 0; i < n; i++){
		c=0;
		for(k = sub*j; k < (j+1)*sub; k++)
		{	
			uj[i*sub+c] = ds[i*d+k];
			c++;
		}
	}																		
	return uj;
}

/**
 * metodo che prende un sottogruppo (sub dimensionale) del data set
 * j serve per prendere il j-esimo gruppetto di dimensioni j=2 equivale ad U2
 * presupponendo ds salvato per colonna*/
MATRIX Uj_x(MATRIX qs, int j,int m,int n,int d){
	int i,k,c;
	int sub=d/m;
	MATRIX uj = alloc_matrix(n,sub);
	c=0;
	for(k = j*sub; k < (j+1)*sub; k++){	
		uj[c] = qs[k];
		c++;
	}
	return uj;
}
/*metodo che calcola il centroide più vicino al parametro x*/
int centXA(float * centroids, float * x, int k, int d){	
	float dis2=0;
	int park2 = 0;
	float tmp2 = 0;
	int k2=k;
	int d2=d;
	cent_XA(centroids,x,k2,d2,&tmp2,&park2,&dis2);
	return park2;
}
/*metodo che calcola il centroide più vicino al parametro x*/
int centXU(float * centroids, float * x, int k, int d){	
	float dis2;
	int park2 = 0;
	float tmp2 = 0;
	int k2=k;
	int d2=d;
	cent_XU(centroids,x,k2,d2,&tmp2,&park2,&dis2);
	return park2;
}


/*k-means Aligned per la ricerca esaustiva*/
void k_means_colA(MATRIX data, int n, int d, int k, float t, int* labels, MATRIX centroids,int t_min,int t_max) {

	printf("\n--------ALIGNED------------------\n");
	float calc;
	float* min_distance = alloc_matrix(size,1);
	float* distance = alloc_matrix(size,1);
	float offset;
	int iter=0;
	int h, i, j, k_p;
	
	/* size of each cluster */
	float* counts = alloc_matrix(k,1);
	float old_error, error = FLT_MAX;
	MATRIX c = centroids;
	
	/* temp centroids */
	MATRIX c1 = alloc_matrix(k,d);
	
	
	for (j = 0; j < d;j++){
		for (i = 0; i < k;i++) {
			c[i*d+j] = data[i+j*n];		
		}
	}
		



	do {
		iter++;
		
		old_error = error, error = 0;
		
		
		
		clearCentroids(counts,c1,k,d);
		
		
		
		
		

		
		printf("identify the closest cluster in %d iteration\n",iter);

		
		
		for (i = 0; i <= n-size; i+=size){  	//per ogni punto del ds
			assignValue(min_distance,&max_f,size);
			


			
			
			for (j = 0; j < k; j++){ // per ogni centroide

				

				colDistance64OptimizedA(data,centroids,distance,i,j,d,n);



				distanceControl64(distance,min_distance,labels,j,i);
			




				
				
			}
		
			
			
		
			
			

			
			
	
			for (j = 0; j < d; j++){
				
				c1[labels[i]*d+j] += data[i+j*n];
				c1[labels[i+1]*d+j] += data[i+1+j*n];
				c1[labels[i+2]*d+j] += data[i+2+j*n];
				c1[labels[i+3]*d+j] += data[i+3+j*n];
				c1[labels[i+4]*d+j] += data[i+4+j*n];
				c1[labels[i+5]*d+j] += data[i+5+j*n];
				c1[labels[i+6]*d+j] += data[i+6+j*n];
				c1[labels[i+7]*d+j] += data[i+7+j*n];
				c1[labels[i+8]*d+j] += data[i+8+j*n];
				c1[labels[i+9]*d+j] += data[i+9+j*n];
				c1[labels[i+10]*d+j] += data[i+10+j*n];
				c1[labels[i+11]*d+j] += data[i+11+j*n];
				c1[labels[i+12]*d+j] += data[i+12+j*n];
				c1[labels[i+13]*d+j] += data[i+13+j*n];
				c1[labels[i+14]*d+j] += data[i+14+j*n];
				c1[labels[i+15]*d+j] += data[i+15+j*n];
				c1[labels[i+16]*d+j] += data[i+16+j*n];
				c1[labels[i+17]*d+j] += data[i+17+j*n];
				c1[labels[i+18]*d+j] += data[i+18+j*n];
				c1[labels[i+19]*d+j] += data[i+19+j*n];
				c1[labels[i+20]*d+j] += data[i+20+j*n];
				c1[labels[i+21]*d+j] += data[i+21+j*n];
				c1[labels[i+22]*d+j] += data[i+22+j*n];
				c1[labels[i+23]*d+j] += data[i+23+j*n];
				c1[labels[i+24]*d+j] += data[i+24+j*n];
				c1[labels[i+25]*d+j] += data[i+25+j*n];
				c1[labels[i+26]*d+j] += data[i+26+j*n];
				c1[labels[i+27]*d+j] += data[i+27+j*n];
				c1[labels[i+28]*d+j] += data[i+28+j*n];
				c1[labels[i+29]*d+j] += data[i+29+j*n];
				c1[labels[i+30]*d+j] += data[i+30+j*n];
				c1[labels[i+31]*d+j] += data[i+31+j*n];
			}
		


			counts[labels[i]]++;
			counts[labels[i+1]]++;
			counts[labels[i+2]]++;
			counts[labels[i+3]]++;
			counts[labels[i+4]]++;
			counts[labels[i+5]]++;
			counts[labels[i+6]]++;
			counts[labels[i+7]]++;
			counts[labels[i+8]]++;
			counts[labels[i+9]]++;
			counts[labels[i+10]]++;
			counts[labels[i+11]]++;
			counts[labels[i+12]]++;
			counts[labels[i+13]]++;
			counts[labels[i+14]]++;
			counts[labels[i+15]]++;
			counts[labels[i+16]]++;
			counts[labels[i+17]]++;
			counts[labels[i+18]]++;
			counts[labels[i+19]]++;
			counts[labels[i+20]]++;
			counts[labels[i+21]]++;
			counts[labels[i+22]]++;
			counts[labels[i+23]]++;
			counts[labels[i+24]]++;
			counts[labels[i+25]]++;
			counts[labels[i+26]]++;
			counts[labels[i+27]]++;
			counts[labels[i+28]]++;
			counts[labels[i+29]]++;
			counts[labels[i+30]]++;
			counts[labels[i+31]]++;


			/* update standard error */

			error += min_distance[0];
			error += min_distance[1];
			error += min_distance[2];
			error += min_distance[3];				
			error += min_distance[4];
			error += min_distance[5];
			error += min_distance[6];
			error += min_distance[7];
			error += min_distance[8];
			error += min_distance[9];
			error += min_distance[10];
			error += min_distance[11];
			error += min_distance[12];
			error += min_distance[13];
			error += min_distance[14];
			error += min_distance[15];
			error += min_distance[16];
			error += min_distance[17];
			error += min_distance[18];
			error += min_distance[19];				
			error += min_distance[20];
			error += min_distance[21];
			error += min_distance[22];
			error += min_distance[23];
			error += min_distance[24];
			error += min_distance[25];
			error += min_distance[26];
			error += min_distance[27];
			error += min_distance[28];
			error += min_distance[29];
			error += min_distance[30];
			error += min_distance[31];

		}
		for (; i <= n-p; i+=p){  	//per ogni punto del ds
			assignValue(min_distance,&max_f,p);


			
			
			for (j = 0; j < k; j++){ // per ogni centroide

				
				
				
				colDistance64A(data,centroids,distance,i,j,d,n);

				

				
				
				

				
				
				
				distanceControl64Sing(distance,min_distance,labels,j,i);
				



				
				
			}
		
			
			
		
			
			

			
			
			for (j = 0; j < d; j++){
				
				c1[labels[i]*d+j] += data[i+j*n];
				c1[labels[i+1]*d+j] += data[i+1+j*n];
				c1[labels[i+2]*d+j] += data[i+2+j*n];
				c1[labels[i+3]*d+j] += data[i+3+j*n];
				c1[labels[i+4]*d+j] += data[i+4+j*n];
				c1[labels[i+5]*d+j] += data[i+5+j*n];
				c1[labels[i+6]*d+j] += data[i+6+j*n];
				c1[labels[i+7]*d+j] += data[i+7+j*n];
			}
		


			counts[labels[i]]++;
			counts[labels[i+1]]++;
			counts[labels[i+2]]++;
			counts[labels[i+3]]++;
			counts[labels[i+4]]++;
			counts[labels[i+5]]++;
			counts[labels[i+6]]++;
			counts[labels[i+7]]++;

			/* update standard error */

			

			error += min_distance[0];
			error += min_distance[1];
			error += min_distance[2];
			error += min_distance[3];	
			error += min_distance[4];
			error += min_distance[5];
			error += min_distance[6];
			error += min_distance[7];				


		
		}
	
		
		updateCentroid(c,c1,counts,k,d);
		
		
		
		
	if(error == old_error)
		calc=0;
	else if(error > old_error){
		calc = (fabs(error-old_error)/error);
		
	}else{
		calc = (fabs(error-old_error)/old_error);
		
	}

	
	}while (!(t_min <= iter && ((t_max < iter) || calc <= t)));

	dealloc_matrix(counts);

	dealloc_matrix(distance);

	dealloc_matrix(min_distance);

}//k_means

/*k-means Aligned per la ricerca non-esaustiva*/
void NE_k_means_colA(MATRIX data, int n, int d, int k, float t, int* labels, MATRIX centroids,int t_min,int t_max,int nr) {

	printf("\n--------ALIGNED------------------\n");

	float calc;
	

	float* min_distance = alloc_matrix(size,1);
	float* distance = alloc_matrix(size,1);
	float offset;
	int iter=0;
	int h, i, j, k_p; 
	
	float* counts = alloc_matrix(k,1);
	float old_error, error = FLT_MAX;
	MATRIX c = centroids;
	
	/* temp centroids */
	MATRIX c1 = alloc_matrix(k,d);
	
	
	for (j = 0; j < d;j++){
		for (i = 0; i < k;i++) {
			c[i*d+j] = data[i+j*n];		
		}
	}
		



	do {
		iter++;
		/* save error from last step */
		old_error = error, error = 0;
		
		
		
		clearCentroids(counts,c1,k,d);
		
		
		
		
		

		
		printf("identify the closest cluster in %d iteration\n",iter);

		
		
		for (i = 0; i <= nr-size; i+=size){  	//per ogni punto del ds
			assignValue(min_distance,&max_f,size);
			


			
			
			for (j = 0; j < k; j++){ // per ogni centroide

				

				colDistance64OptimizedA(data,centroids,distance,i,j,d,n);



				

				
				
				

				
				
				
				distanceControl64(distance,min_distance,labels,j,i);
			




				
				
			}
		
			
			
		
			
			

			
			
	
			for (j = 0; j < d; j++){
				
				c1[labels[i]*d+j] += data[i+j*n];
				c1[labels[i+1]*d+j] += data[i+1+j*n];
				c1[labels[i+2]*d+j] += data[i+2+j*n];
				c1[labels[i+3]*d+j] += data[i+3+j*n];
				c1[labels[i+4]*d+j] += data[i+4+j*n];
				c1[labels[i+5]*d+j] += data[i+5+j*n];
				c1[labels[i+6]*d+j] += data[i+6+j*n];
				c1[labels[i+7]*d+j] += data[i+7+j*n];
				c1[labels[i+8]*d+j] += data[i+8+j*n];
				c1[labels[i+9]*d+j] += data[i+9+j*n];
				c1[labels[i+10]*d+j] += data[i+10+j*n];
				c1[labels[i+11]*d+j] += data[i+11+j*n];
				c1[labels[i+12]*d+j] += data[i+12+j*n];
				c1[labels[i+13]*d+j] += data[i+13+j*n];
				c1[labels[i+14]*d+j] += data[i+14+j*n];
				c1[labels[i+15]*d+j] += data[i+15+j*n];
				c1[labels[i+16]*d+j] += data[i+16+j*n];
				c1[labels[i+17]*d+j] += data[i+17+j*n];
				c1[labels[i+18]*d+j] += data[i+18+j*n];
				c1[labels[i+19]*d+j] += data[i+19+j*n];
				c1[labels[i+20]*d+j] += data[i+20+j*n];
				c1[labels[i+21]*d+j] += data[i+21+j*n];
				c1[labels[i+22]*d+j] += data[i+22+j*n];
				c1[labels[i+23]*d+j] += data[i+23+j*n];
				c1[labels[i+24]*d+j] += data[i+24+j*n];
				c1[labels[i+25]*d+j] += data[i+25+j*n];
				c1[labels[i+26]*d+j] += data[i+26+j*n];
				c1[labels[i+27]*d+j] += data[i+27+j*n];
				c1[labels[i+28]*d+j] += data[i+28+j*n];
				c1[labels[i+29]*d+j] += data[i+29+j*n];
				c1[labels[i+30]*d+j] += data[i+30+j*n];
				c1[labels[i+31]*d+j] += data[i+31+j*n];
			}
		


			counts[labels[i]]++;
			counts[labels[i+1]]++;
			counts[labels[i+2]]++;
			counts[labels[i+3]]++;
			counts[labels[i+4]]++;
			counts[labels[i+5]]++;
			counts[labels[i+6]]++;
			counts[labels[i+7]]++;
			counts[labels[i+8]]++;
			counts[labels[i+9]]++;
			counts[labels[i+10]]++;
			counts[labels[i+11]]++;
			counts[labels[i+12]]++;
			counts[labels[i+13]]++;
			counts[labels[i+14]]++;
			counts[labels[i+15]]++;
			counts[labels[i+16]]++;
			counts[labels[i+17]]++;
			counts[labels[i+18]]++;
			counts[labels[i+19]]++;
			counts[labels[i+20]]++;
			counts[labels[i+21]]++;
			counts[labels[i+22]]++;
			counts[labels[i+23]]++;
			counts[labels[i+24]]++;
			counts[labels[i+25]]++;
			counts[labels[i+26]]++;
			counts[labels[i+27]]++;
			counts[labels[i+28]]++;
			counts[labels[i+29]]++;
			counts[labels[i+30]]++;
			counts[labels[i+31]]++;
			/* update standard error */

			

			error += min_distance[0];
			error += min_distance[1];
			error += min_distance[2];
			error += min_distance[3];				
			error += min_distance[4];
			error += min_distance[5];
			error += min_distance[6];
			error += min_distance[7];
			error += min_distance[8];
			error += min_distance[9];
			error += min_distance[10];
			error += min_distance[11];
			error += min_distance[12];
			error += min_distance[13];
			error += min_distance[14];
			error += min_distance[15];
			error += min_distance[16];
			error += min_distance[17];
			error += min_distance[18];
			error += min_distance[19];				
			error += min_distance[20];
			error += min_distance[21];
			error += min_distance[22];
			error += min_distance[23];
			error += min_distance[24];
			error += min_distance[25];
			error += min_distance[26];
			error += min_distance[27];
			error += min_distance[28];
			error += min_distance[29];
			error += min_distance[30];
			error += min_distance[31];

		}
		for (; i <= nr-p; i+=p){  	//per ogni punto del ds
			assignValue(min_distance,&max_f,p);


			
			
			for (j = 0; j < k; j++){ // per ogni centroide

				
				
				
				colDistance64A(data,centroids,distance,i,j,d,n);

				

				
				
				

				
				
				
				distanceControl64Sing(distance,min_distance,labels,j,i);
				



				
				
			}
		
			
			
		
			
			

			
			
			for (j = 0; j < d; j++){
				
				c1[labels[i]*d+j] += data[i+j*n];
				c1[labels[i+1]*d+j] += data[i+1+j*n];
				c1[labels[i+2]*d+j] += data[i+2+j*n];
				c1[labels[i+3]*d+j] += data[i+3+j*n];
				c1[labels[i+4]*d+j] += data[i+4+j*n];
				c1[labels[i+5]*d+j] += data[i+5+j*n];
				c1[labels[i+6]*d+j] += data[i+6+j*n];
				c1[labels[i+7]*d+j] += data[i+7+j*n];
			}
		


			counts[labels[i]]++;
			counts[labels[i+1]]++;
			counts[labels[i+2]]++;
			counts[labels[i+3]]++;
			counts[labels[i+4]]++;
			counts[labels[i+5]]++;
			counts[labels[i+6]]++;
			counts[labels[i+7]]++;

			/* update standard error */

			

			error += min_distance[0];
			error += min_distance[1];
			error += min_distance[2];
			error += min_distance[3];	
			error += min_distance[4];
			error += min_distance[5];
			error += min_distance[6];
			error += min_distance[7];				


		
		}
	
		for (; i < nr; i++){  	//per ogni punto del ds
			assignValue(min_distance,&max_f,p);


			
			
			for (j = 0; j < k; j++){ // per ogni centroide

				
			
				
				colDistance64Sing(data,centroids,distance,i,j,d,n);
				
				
				

				if (distance[0] < min_distance[0]) {
					labels[i] = j;
					min_distance[0] = distance[0];
				}
				
				

				
				



				
				
			}
		
			
			
		
			
			

			
			
			for (j = 0; j < d; j++){
				
				c1[labels[i]*d+j] += data[i+j*n];

			}
		


			counts[labels[i]]++;

			/* update standard error */

			

			error += min_distance[0];
				


		
		}

		updateCentroid(c,c1,counts,k,d);
		
		
		
		
		
	if(error == old_error)
		calc=0;
	else if(error > old_error){
		calc = (fabs(error-old_error)/error);
		
	}else{
		calc = (fabs(error-old_error)/old_error);

	}

	
	}while (!(t_min <= iter && ((t_max < iter) || calc <= t)));


		
	for (; i <= n-size; i+=size){  	//per ogni punto del ds
		assignValue(min_distance,&max_f,size);
		


		
		
		for (j = 0; j < k; j++){ // per ogni centroide

			

			colDistance64OptimizedA(data,centroids,distance,i,j,d,n);



			
			
			distanceControl64(distance,min_distance,labels,j,i);
		



		

			
			
		}
	
		
		
	
		

		


	}
	for (; i <= n-p; i+=p){  	//per ogni punto del ds
		assignValue(min_distance,&max_f,p);


		
		
		for (j = 0; j < k; j++){ // per ogni centroide

			
			
			
			colDistance64A(data,centroids,distance,i,j,d,n);

			
			
			

		
			
			
			distanceControl64Sing(distance,min_distance,labels,j,i);
			


		

			
			
		}
	

	
	}
	for (; i < n; i++){  	//per ogni punto del ds
		
		assignValue(min_distance,&max_f,p);


		
		
		for (j = 0; j < k; j++){ // per ogni centroide

			
			colDistance64Sing(data,centroids,distance,i,j,d,n);
			
			if (distance[0] < min_distance[0]) {
				labels[i] = j;
				min_distance[0] = distance[0];
			}
			
			

		

		

			
			
		}

	}
	
	dealloc_matrix(counts);

	dealloc_matrix(distance);

	dealloc_matrix(min_distance);

}//k_means

/*k-means Unaligned per la ricerca esaustiva*/
void k_means_colU(MATRIX data, int n, int d, int k, float t, int* labels, MATRIX centroids,int t_min,int t_max) {

	printf("\n--------UNALIGNED------------------\n");


	
	float calc;
	float* min_distance = alloc_matrix(size,1);
	float* distance = alloc_matrix(size,1);
	float offset;
	int iter=0;
	int h, i, j, k_p; 
	
	float* counts = alloc_matrix(k,1);
	float old_error, error = FLT_MAX;
	MATRIX c = centroids;
	
	MATRIX c1 = alloc_matrix(k,d);
	
	

	for (j = 0; j < d;j++){
		
		for (i = 0; i < k ;i++) {
			c[i*d+j] = data[i+j*n];		
		}
	}
	

	do {
		iter++;
		/* save error from last step */
		old_error = error, error = 0;
		
		
		
		clearCentroids(counts,c1,k,d);
		
		

		
		

		
		printf("identify the closest cluster in %d iteration\n",iter);

		
		
		for (i = 0; i <= n-size; i+=size){  	//per ogni punto del ds
			assignValue(min_distance,&max_f,size);


			
			
			for (j = 0; j < k; j++){ // per ogni centroide

				
				colDistance64OptimizedU(data,centroids,distance,i,j,d,n);



				
				

				
				
				
				distanceControl64(distance,min_distance,labels,j,i);
			




				
				
			}
		
			
			
		
			
			

			
			
	
			for (j = 0; j < d; j++){
				
				c1[labels[i]*d+j] += data[i+j*n];
				c1[labels[i+1]*d+j] += data[i+1+j*n];
				c1[labels[i+2]*d+j] += data[i+2+j*n];
				c1[labels[i+3]*d+j] += data[i+3+j*n];
				c1[labels[i+4]*d+j] += data[i+4+j*n];
				c1[labels[i+5]*d+j] += data[i+5+j*n];
				c1[labels[i+6]*d+j] += data[i+6+j*n];
				c1[labels[i+7]*d+j] += data[i+7+j*n];
				c1[labels[i+8]*d+j] += data[i+8+j*n];
				c1[labels[i+9]*d+j] += data[i+9+j*n];
				c1[labels[i+10]*d+j] += data[i+10+j*n];
				c1[labels[i+11]*d+j] += data[i+11+j*n];
				c1[labels[i+12]*d+j] += data[i+12+j*n];
				c1[labels[i+13]*d+j] += data[i+13+j*n];
				c1[labels[i+14]*d+j] += data[i+14+j*n];
				c1[labels[i+15]*d+j] += data[i+15+j*n];
				c1[labels[i+16]*d+j] += data[i+16+j*n];
				c1[labels[i+17]*d+j] += data[i+17+j*n];
				c1[labels[i+18]*d+j] += data[i+18+j*n];
				c1[labels[i+19]*d+j] += data[i+19+j*n];
				c1[labels[i+20]*d+j] += data[i+20+j*n];
				c1[labels[i+21]*d+j] += data[i+21+j*n];
				c1[labels[i+22]*d+j] += data[i+22+j*n];
				c1[labels[i+23]*d+j] += data[i+23+j*n];
				c1[labels[i+24]*d+j] += data[i+24+j*n];
				c1[labels[i+25]*d+j] += data[i+25+j*n];
				c1[labels[i+26]*d+j] += data[i+26+j*n];
				c1[labels[i+27]*d+j] += data[i+27+j*n];
				c1[labels[i+28]*d+j] += data[i+28+j*n];
				c1[labels[i+29]*d+j] += data[i+29+j*n];
				c1[labels[i+30]*d+j] += data[i+30+j*n];
				c1[labels[i+31]*d+j] += data[i+31+j*n];
			}
		


			counts[labels[i]]++;
			counts[labels[i+1]]++;
			counts[labels[i+2]]++;
			counts[labels[i+3]]++;
			counts[labels[i+4]]++;
			counts[labels[i+5]]++;
			counts[labels[i+6]]++;
			counts[labels[i+7]]++;
			counts[labels[i+8]]++;
			counts[labels[i+9]]++;
			counts[labels[i+10]]++;
			counts[labels[i+11]]++;
			counts[labels[i+12]]++;
			counts[labels[i+13]]++;
			counts[labels[i+14]]++;
			counts[labels[i+15]]++;
			counts[labels[i+16]]++;
			counts[labels[i+17]]++;
			counts[labels[i+18]]++;
			counts[labels[i+19]]++;
			counts[labels[i+20]]++;
			counts[labels[i+21]]++;
			counts[labels[i+22]]++;
			counts[labels[i+23]]++;
			counts[labels[i+24]]++;
			counts[labels[i+25]]++;
			counts[labels[i+26]]++;
			counts[labels[i+27]]++;
			counts[labels[i+28]]++;
			counts[labels[i+29]]++;
			counts[labels[i+30]]++;
			counts[labels[i+31]]++;
			/* update standard error */

			

			error += min_distance[0];
			error += min_distance[1];
			error += min_distance[2];
			error += min_distance[3];				
			error += min_distance[4];
			error += min_distance[5];
			error += min_distance[6];
			error += min_distance[7];
			error += min_distance[8];
			error += min_distance[9];
			error += min_distance[10];
			error += min_distance[11];
			error += min_distance[12];
			error += min_distance[13];
			error += min_distance[14];
			error += min_distance[15];
			error += min_distance[16];
			error += min_distance[17];
			error += min_distance[18];
			error += min_distance[19];				
			error += min_distance[20];
			error += min_distance[21];
			error += min_distance[22];
			error += min_distance[23];
			error += min_distance[24];
			error += min_distance[25];
			error += min_distance[26];
			error += min_distance[27];
			error += min_distance[28];
			error += min_distance[29];
			error += min_distance[30];
			error += min_distance[31];


		}
		for (; i <= n-p; i+=p){  	//per ogni punto del ds
			assignValue(min_distance,&max_f,p);


			
			
			for (j = 0; j < k; j++){ // per ogni centroide

				
			
				
				colDistance64U(data,centroids,distance,i,j,d,n);

				

				
				
				

				
				
				
				distanceControl64Sing(distance,min_distance,labels,j,i);
				



				
				
			}
		
			
			
		
			
			

			
			
			for (j = 0; j < d; j++){
				
				c1[labels[i]*d+j] += data[i+j*n];
				c1[labels[i+1]*d+j] += data[i+1+j*n];
				c1[labels[i+2]*d+j] += data[i+2+j*n];
				c1[labels[i+3]*d+j] += data[i+3+j*n];
				c1[labels[i+4]*d+j] += data[i+4+j*n];
				c1[labels[i+5]*d+j] += data[i+5+j*n];
				c1[labels[i+6]*d+j] += data[i+6+j*n];
				c1[labels[i+7]*d+j] += data[i+7+j*n];
			}
		


			counts[labels[i]]++;
			counts[labels[i+1]]++;
			counts[labels[i+2]]++;
			counts[labels[i+3]]++;
			counts[labels[i+4]]++;
			counts[labels[i+5]]++;
			counts[labels[i+6]]++;
			counts[labels[i+7]]++;

			/* update standard error */

			

			error += min_distance[0];
			error += min_distance[1];
			error += min_distance[2];
			error += min_distance[3];	
			error += min_distance[4];
			error += min_distance[5];
			error += min_distance[6];
			error += min_distance[7];				


		
		}
		

		for (; i < n; i++){  	//per ogni punto del ds
			assignValue(min_distance,&max_f,p);


			
			
			for (j = 0; j < k; j++){ // per ogni centroide

				
			
				
				colDistance64Sing(data,centroids,distance,i,j,d,n);
				
				
				

				if (distance[0] < min_distance[0]) {
					labels[i] = j;
					min_distance[0] = distance[0];
				}
				
				

				
				



				
				
			}
		
			
			
		
			
			

			
			
			for (j = 0; j < d; j++){
				
				c1[labels[i]*d+j] += data[i+j*n];

			}
		


			counts[labels[i]]++;

			/* update standard error */

			

			error += min_distance[0];
				


		
		}


		
		

		
		
	
		
		updateCentroid(c,c1,counts,k,d);
		
		
		
		
		
	if(error == old_error)
		calc=0;
	else if(error > old_error){
		calc = (fabs(error-old_error)/error);
		
	}else{
		calc = (fabs(error-old_error)/old_error);
		
	}

	
	}while (!(t_min <= iter && ((t_max < iter) || calc <= t)));

	

	dealloc_matrix(counts);

	dealloc_matrix(distance);

	dealloc_matrix(min_distance);

}//k_means
/*k-means Unaligned per la ricerca non-esaustiva*/
void NE_k_means_colU(MATRIX data, int n, int d, int k, float t, int* labels, MATRIX centroids,int t_min,int t_max,int nr) {

	printf("\n--------UNALIGNED------------------\n");


	
	float calc;
	float* min_distance = alloc_matrix(size,1);
	float* distance = alloc_matrix(size,1);
	float offset;
	int iter=0;
	int h, i, j, k_p; 
	
	float* counts = alloc_matrix(k,1);
	float old_error, error = FLT_MAX;
	MATRIX c = centroids;
	
	MATRIX c1 = alloc_matrix(k,d);
	
	

	for (j = 0; j < d;j++){
		
		for (i = 0; i < k ;i++) {
			c[i*d+j] = data[i+j*n];		
		}
	}
	
	do {
		iter++;
		/* save error from last step */
		old_error = error, error = 0;
		
		
		
		clearCentroids(counts,c1,k,d);
		
		

		
		

		
		printf("identify the closest cluster in %d iteration\n",iter);

		
		
		for (i = 0; i <= nr-size; i+=size){  	//per ogni punto del ds
			assignValue(min_distance,&max_f,size);


			
			
			for (j = 0; j < k; j++){ // per ogni centroide

				
				colDistance64OptimizedU(data,centroids,distance,i,j,d,n);



				
				

				
				
				
				distanceControl64(distance,min_distance,labels,j,i);
			




				
				
			}
		
			
			
		
			
			

			
			
	
			for (j = 0; j < d; j++){
				
				c1[labels[i]*d+j] += data[i+j*n];
				c1[labels[i+1]*d+j] += data[i+1+j*n];
				c1[labels[i+2]*d+j] += data[i+2+j*n];
				c1[labels[i+3]*d+j] += data[i+3+j*n];
				c1[labels[i+4]*d+j] += data[i+4+j*n];
				c1[labels[i+5]*d+j] += data[i+5+j*n];
				c1[labels[i+6]*d+j] += data[i+6+j*n];
				c1[labels[i+7]*d+j] += data[i+7+j*n];
				c1[labels[i+8]*d+j] += data[i+8+j*n];
				c1[labels[i+9]*d+j] += data[i+9+j*n];
				c1[labels[i+10]*d+j] += data[i+10+j*n];
				c1[labels[i+11]*d+j] += data[i+11+j*n];
				c1[labels[i+12]*d+j] += data[i+12+j*n];
				c1[labels[i+13]*d+j] += data[i+13+j*n];
				c1[labels[i+14]*d+j] += data[i+14+j*n];
				c1[labels[i+15]*d+j] += data[i+15+j*n];
				c1[labels[i+16]*d+j] += data[i+16+j*n];
				c1[labels[i+17]*d+j] += data[i+17+j*n];
				c1[labels[i+18]*d+j] += data[i+18+j*n];
				c1[labels[i+19]*d+j] += data[i+19+j*n];
				c1[labels[i+20]*d+j] += data[i+20+j*n];
				c1[labels[i+21]*d+j] += data[i+21+j*n];
				c1[labels[i+22]*d+j] += data[i+22+j*n];
				c1[labels[i+23]*d+j] += data[i+23+j*n];
				c1[labels[i+24]*d+j] += data[i+24+j*n];
				c1[labels[i+25]*d+j] += data[i+25+j*n];
				c1[labels[i+26]*d+j] += data[i+26+j*n];
				c1[labels[i+27]*d+j] += data[i+27+j*n];
				c1[labels[i+28]*d+j] += data[i+28+j*n];
				c1[labels[i+29]*d+j] += data[i+29+j*n];
				c1[labels[i+30]*d+j] += data[i+30+j*n];
				c1[labels[i+31]*d+j] += data[i+31+j*n];
			}
		


			counts[labels[i]]++;
			counts[labels[i+1]]++;
			counts[labels[i+2]]++;
			counts[labels[i+3]]++;
			counts[labels[i+4]]++;
			counts[labels[i+5]]++;
			counts[labels[i+6]]++;
			counts[labels[i+7]]++;
			counts[labels[i+8]]++;
			counts[labels[i+9]]++;
			counts[labels[i+10]]++;
			counts[labels[i+11]]++;
			counts[labels[i+12]]++;
			counts[labels[i+13]]++;
			counts[labels[i+14]]++;
			counts[labels[i+15]]++;
			counts[labels[i+16]]++;
			counts[labels[i+17]]++;
			counts[labels[i+18]]++;
			counts[labels[i+19]]++;
			counts[labels[i+20]]++;
			counts[labels[i+21]]++;
			counts[labels[i+22]]++;
			counts[labels[i+23]]++;
			counts[labels[i+24]]++;
			counts[labels[i+25]]++;
			counts[labels[i+26]]++;
			counts[labels[i+27]]++;
			counts[labels[i+28]]++;
			counts[labels[i+29]]++;
			counts[labels[i+30]]++;
			counts[labels[i+31]]++;
			/* update standard error */

			

			error += min_distance[0];
			error += min_distance[1];
			error += min_distance[2];
			error += min_distance[3];				
			error += min_distance[4];
			error += min_distance[5];
			error += min_distance[6];
			error += min_distance[7];
			error += min_distance[8];
			error += min_distance[9];
			error += min_distance[10];
			error += min_distance[11];
			error += min_distance[12];
			error += min_distance[13];
			error += min_distance[14];
			error += min_distance[15];
			error += min_distance[16];
			error += min_distance[17];
			error += min_distance[18];
			error += min_distance[19];				
			error += min_distance[20];
			error += min_distance[21];
			error += min_distance[22];
			error += min_distance[23];
			error += min_distance[24];
			error += min_distance[25];
			error += min_distance[26];
			error += min_distance[27];
			error += min_distance[28];
			error += min_distance[29];
			error += min_distance[30];
			error += min_distance[31];


		}
		for (; i <= nr-p; i+=p){  	//per ogni punto del ds
			assignValue(min_distance,&max_f,p);


			
			
			for (j = 0; j < k; j++){ // per ogni centroide

				
			
				
				colDistance64U(data,centroids,distance,i,j,d,n);

				

				
				
				

				
				
				
				distanceControl64Sing(distance,min_distance,labels,j,i);
				



				
				
			}
		
			
			
		
			
			

			
			
			for (j = 0; j < d; j++){
				
				c1[labels[i]*d+j] += data[i+j*n];
				c1[labels[i+1]*d+j] += data[i+1+j*n];
				c1[labels[i+2]*d+j] += data[i+2+j*n];
				c1[labels[i+3]*d+j] += data[i+3+j*n];
				c1[labels[i+4]*d+j] += data[i+4+j*n];
				c1[labels[i+5]*d+j] += data[i+5+j*n];
				c1[labels[i+6]*d+j] += data[i+6+j*n];
				c1[labels[i+7]*d+j] += data[i+7+j*n];
			}
		


			counts[labels[i]]++;
			counts[labels[i+1]]++;
			counts[labels[i+2]]++;
			counts[labels[i+3]]++;
			counts[labels[i+4]]++;
			counts[labels[i+5]]++;
			counts[labels[i+6]]++;
			counts[labels[i+7]]++;

			/* update standard error */

			

			error += min_distance[0];
			error += min_distance[1];
			error += min_distance[2];
			error += min_distance[3];	
			error += min_distance[4];
			error += min_distance[5];
			error += min_distance[6];
			error += min_distance[7];				


		
		}
		

		for (; i < nr; i++){  	//per ogni punto del ds
			assignValue(min_distance,&max_f,p);


			
			
			for (j = 0; j < k; j++){ // per ogni centroide

				
			
				
				colDistance64Sing(data,centroids,distance,i,j,d,n);
				
				
				

				if (distance[0] < min_distance[0]) {
					labels[i] = j;
					min_distance[0] = distance[0];
				}
				
				

				
				



				
				
			}
		
			
			
		
			
			

			
			
			for (j = 0; j < d; j++){
				
				c1[labels[i]*d+j] += data[i+j*n];

			}
		


			counts[labels[i]]++;

			/* update standard error */

			

			error += min_distance[0];
				


		
		}


		
		

		
		
	
		
		updateCentroid(c,c1,counts,k,d);
		
		
		
		
		
	if(error == old_error)
		calc=0;
	else if(error > old_error){
		calc = (fabs(error-old_error)/error);
		
	}else{
		calc = (fabs(error-old_error)/old_error);
		
	}

	
	}while (!(t_min <= iter && ((t_max < iter) || calc <= t)));


	for (; i <= n-size; i+=size){  	//per ogni punto del ds
		assignValue(min_distance,&max_f,size);


		
		
		for (j = 0; j < k; j++){ // per ogni centroide

			
			colDistance64OptimizedU(data,centroids,distance,i,j,d,n);




			
			
			

		
			
			
			distanceControl64(distance,min_distance,labels,j,i);
		



		

			
			
		}
	
		
		


	}
	for (; i <= n-p; i+=p){  	//per ogni punto del ds
		assignValue(min_distance,&max_f,p);


		
		
		for (j = 0; j < k; j++){ // per ogni centroide

			
			
			
			colDistance64U(data,centroids,distance,i,j,d,n);

			

			
			
			

		
			
			
			distanceControl64Sing(distance,min_distance,labels,j,i);
			


		

			
			
		}
	
		
		
	
		

				


	
	}
	

	for (; i < n; i++){  	//per ogni punto del ds
		assignValue(min_distance,&max_f,p);


		
		
		for (j = 0; j < k; j++){ // per ogni centroide

			
			
			
			colDistance64Sing(data,centroids,distance,i,j,d,n);
			
			
			if (distance[0] < min_distance[0]) {
				labels[i] = j;
				min_distance[0] = distance[0];
			}
			
			

		
			
			
		}
	
		
		
	
		

	}


	dealloc_matrix(counts);

	dealloc_matrix(distance);

	dealloc_matrix(min_distance);

}//k_means



/*metodo che calcola i residui di x presupponendo il dataset per riga, utilizzato nel calcolo dei residui del query-set*/
MATRIX residuals_row(MATRIX ds,MATRIX centroids,int* label, int n,int d){
	MATRIX results = alloc_matrix(n,d);
	int i,j;
	for (i=0;i< n;i++){
		for(j=0;j<d;j++){
			results[i*d+j]=ds[i*d+j]-centroids[label[i]*d+j];
		}
	}
	
	return results;
}
/*metodo che calcola i residui di x presupponendo il dataset per colonna, utilizzato nel calcolo dei residui del dataset*/
MATRIX residuals_col(MATRIX ds,MATRIX centroids,int* label, int n,int d){
	MATRIX results = alloc_matrix(n,d);
	int i,j;
	for (i=0;i< n;i++){
		for(j=0;j<d;j++){
			results[i+n*j]=ds[i+j*n]-centroids[label[i]*d+j];
		}
	}
	return results;
}

/*implementatione del Producti quantizzation in cui vengono richiamate opportunamente più k_means su ogni sottogruppetto,utilizzato nella variante Aligned con ricerca esaustiva*/
int* productQuantA(MATRIX ds,int n,int d,int m,int k,float* centroids,float eps,int t_min,int t_max){
	int j;
	int sub=d/m;
	int* result = alloc_vector(m*n);
	MATRIX tmp;
	#pragma omp parallel for	
	for( j = 0; j < m; j++)
	{
		if(m!=1){
			tmp = Uj_col(ds,j,m,n,d);
			
			k_means_colA(tmp,n,sub,k,eps,&result[j*n],&centroids[j*sub*k],t_min,t_max);
			
		}else
		{
			k_means_colA(ds,n,sub,k,eps,&result[j*n],&centroids[j*sub*k],t_min,t_max);
		}
	}
	dealloc_matrix(tmp); 
	return result;
}
/*implementatione del Producti quantizzation in cui vengono richiamate opportunamente più k_means su ogni sottogruppetto,utilizzato nella variante Aligned con ricerca non-esaustiva*/
int* NE_productQuantA(MATRIX ds,int n,int d,int m,int k,float* centroids,float eps,int t_min,int t_max,int nr){
	int j;
	int sub=d/m;
	int* result = alloc_vector(m*n);
	MATRIX tmp;
	#pragma omp parallel for	
	for( j = 0; j < m; j++)
	{
		if(m!=1){
			tmp = Uj_col(ds,j,m,n,d);
			
			NE_k_means_colA(tmp,n,sub,k,eps,&result[j*n],&centroids[j*sub*k],t_min,t_max,nr);
			
		}else
		{
			NE_k_means_colA(ds,n,sub,k,eps,&result[j*n],&centroids[j*sub*k],t_min,t_max,nr);
		}
	}
	dealloc_matrix(tmp); 
	return result;
}
/*implementatione del Producti quantizzation in cui vengono richiamate opportunamente più k_means su ogni sottogruppetto,utilizzato nella variante Unaligned con ricerca esaustiva*/
int* productQuantU(MATRIX ds,int n,int d,int m,int k,float* centroids,float eps,int t_min,int t_max){
	int j;
	int sub=d/m;
	int* result = alloc_vector(m*n);
	MATRIX tmp;
	#pragma omp parallel for	
	for( j = 0; j < m; j++)
	{
		if(m!=1){
			tmp = Uj_col(ds,j,m,n,d);
			
			k_means_colU(tmp,n,sub,k,eps,&result[j*n],&centroids[j*sub*k],t_min,t_max);
			
		}else
		{
			k_means_colU(ds,n,sub,k,eps,&result[j*n],&centroids[j*sub*k],t_min,t_max);
		}
	}
	dealloc_matrix(tmp); 
	return result;
}
/*implementatione del Producti quantizzation in cui vengono richiamate opportunamente più k_means su ogni sottogruppetto,utilizzato nella variante Unaligned con ricerca non-esaustiva*/
int* NE_productQuantU(MATRIX ds,int n,int d,int m,int k,float* centroids,float eps,int t_min,int t_max,int nr){
	int j;
	int sub=d/m;
	int* result = alloc_vector(m*n);
	MATRIX tmp;
	#pragma omp parallel for
	for( j = 0; j < m; j++)
	{
		if(m!=1){
			tmp = Uj_col(ds,j,m,n,d);
			
			NE_k_means_colU(tmp,n,sub,k,eps,&result[j*n],&centroids[j*sub*k],t_min,t_max,nr);
			
		}else
		{
			NE_k_means_colU(ds,n,sub,k,eps,&result[j*n],&centroids[j*sub*k],t_min,t_max,nr);
		}
	}
	dealloc_matrix(tmp); 
	return result;
}
/*metodo che gestisce la struttura max_heap, struttura tenuta non ordinata e che esporta una variabile che mantiene il massimo corrente all'interno della struttura*/
void max_heap(int* index,float* result_dist,int y,float tmp,float max,int dim,bool full, float* result){
	float new_max;
	if (full==true || c_max_heap == dim) { 
		bool trovato = false;
		new_max=tmp;
		for(int j = 0; j < dim; j++)
		{
			
			if(!trovato && result_dist[j]==max){
				result_dist[j]=tmp;
				index[j]=y;
				trovato=true;
			}else if (result_dist[j] > new_max)
			{
				new_max = result_dist[j];
			}	
		}
		result[0] = new_max;
	}
	else
	{
		index[c_max_heap]=y;
		result_dist[c_max_heap]=tmp;
		if (tmp > pre_max_heap) {
			pre_max_heap=tmp;
		}
		c_max_heap++;
		
		if (c_max_heap==dim) {
			result[0] = pre_max_heap;
		}else
		{
			result[0] = max;
		}
	}

}


/*metodo che calcola i w centroidi più vicini ad x, utilizzato nella ricerca non esaustiva. Variante Aligned*/
int * w_near_centroidsA(MATRIX x,MATRIX centroids,int n,int d,int w){
	int i,j;
	int * result_w=alloc_vector(w);
	float * result_dist=alloc_matrix(w,1);
	float tmp=0;
	float max=0;
	int new_i;
	
	for(i = 0; i < w; i++)
	{	
		tmp = 0;
		dist64A(x, &centroids[i*d], &tmp, d);
		
		result_w[i]=i;
		result_dist[i]=tmp;
	
		if (max < tmp) {
			max = tmp;
		}
	}

	

	for(i=w;i<n;i++){
	
		tmp = 0;
		dist64A(x, &centroids[i*d], &tmp, d);
	
		
		if(tmp < max){

			max_heap(result_w,result_dist,i,tmp,max,w,true, &max);

		}

	}
	dealloc_matrix(result_dist);
	return result_w;


}
/*metodo che calcola i w centroidi più vicini ad x, utilizzato nella ricerca non esaustiva. Variante Unaligned*/
int * w_near_centroidsU(MATRIX x,MATRIX centroids,int n,int d,int w){
	int i,j;
	int * result_w=alloc_vector(w);
	float * result_dist=alloc_matrix(w,1);
	float tmp=0;
	float max=0;
	int new_i;

	for(i = 0; i < w; i++)
	{	
		tmp = 0;
		dist64U(x, &centroids[i*d], &tmp, d);
	
		result_w[i]=i;
		result_dist[i]=tmp;
	
		if (max < tmp) {
			max = tmp;
		}
	}

	

	for(i=w;i<n;i++){
	
		tmp = 0;
		dist64U(x, &centroids[i*d], &tmp, d);
	
		
		if(tmp < max){
	
			max_heap(result_w,result_dist,i,tmp,max,w,true, &max);
	
			
		}

	}
	dealloc_matrix(result_dist);
	return result_w;


}
/*metodo per il calcolo dei residui di X*/
MATRIX residuals_x(MATRIX x,MATRIX centroids,int* label, int n,int d){
	MATRIX results = alloc_matrix(n,d);
	int i,j;
	for (i=0;i< n;i++){
		for(j=0;j<d;j++){
			results[i*d+j]=x[j]-centroids[label[i]*d+j];
		}
	}
	return results;
}


/*metodo per il precalcolo delle distanze come descritto in ADC. Variante Aligned*/
float* pre_adcA(MATRIX x, float* centroids,int d,int m, int k, int sub ){
	float* result= alloc_matrix(m,k);
	int i,j;
	float distance;
	MATRIX uj_x;
	for(j=0; j<m; j++){
		uj_x = Uj_x( x, j, m, 1, d);
	
		for(i = 0; i < k; i++){
			distance = 0;
			rowDistance64AdcA(centroids,uj_x,&distance,i,j,k,sub);
				
			result[j*k+i] = distance;
		

		}
		dealloc_matrix(uj_x);
	}
	
	return result;
}
/*metodo per il precalcolo delle distanze come descritto in ADC. Variante Unaligned*/
float* pre_adcU(MATRIX x, float* centroids,int d,int m, int k, int sub ){
	float* result= alloc_matrix(m,k);
	int i,j;
	float distance,distance2;
	MATRIX uj_x;
	for(j=0; j<m; j++){
		uj_x = Uj_x( x, j, m, 1, d);
	
		for(i = 0; i < k; i++){
			distance = 0;
			rowDistance64AdcU(centroids,uj_x,&distance,i,j,k,sub);
				
			result[j*k+i] = distance;
		

		}
		dealloc_matrix(uj_x);
	}
	return result;
}
/*metodo per il precalcolo delle distanze come descritto in SDC. Variante Aligned*/
float* pre_sdcA(float* centroids,int sub,int m, int k){
	int k_2 = k*k;
	float* result= alloc_matrix(m,k_2);
	int i,j,c,j_d,j_k,i_k;
	float distance;
	for(j=0; j<m; j++){
		j_k = j*k_2;
		for(i = 0; i < k; i++){
			i_k = i*k;
			for(j_d = i+1; j_d < k;j_d++){
				distance = 0;
				rowDistance64SdcA(centroids,&distance,i,j,j_d,k,sub);
				result[j_k+i_k+j_d]=distance;
				result[j_k+j_d*k+i]=distance;
				result[j_k+i_k+i]=0;			
			}
		}
	}
	return result;
}
/*metodo per il precalcolo delle distanze come descritto in SDC. Variante Unaligned*/
float* pre_sdcU(float* centroids,int sub,int m, int k ){
	int k_2 = k*k;
	float* result= alloc_matrix(m,k_2);
	int i,j,c,j_d,j_k,i_k;
	float distance;
	for(j=0; j<m; j++){
		j_k = j*k_2;
		for(i = 0; i < k; i++){
			i_k = i*k;
			for(j_d = i+1; j_d < k;j_d++){
				distance = 0;
				rowDistance64SdcU(centroids,&distance,i,j,j_d,k,sub);
				result[j_k+i_k+j_d]=distance;
				result[j_k+j_d*k+i]=distance;
				result[j_k+i_k+i]=0;			
			}
		}
	}
	return result;
}



void pqnn_index(params* input) {
	int i,j;


	if(input->exaustive == 0 && nmod4==true){


	

		
		Cc= alloc_matrix(input->kc,input->d);
		Cc_index = alloc_vector(input->n);
		NE_k_means_colA(input->ds,input->n,input->d,input->kc,input->eps,Cc_index,Cc,input->tmin,input->tmax,input->nr);
		
		
		printf("Calcolo dei residui\n");
		
		MATRIX res= residuals_col(input->ds,Cc,Cc_index,input->n,input->d);

		
		printf("Quantizzazione dei residui\n");
	

		Cp = alloc_matrix(input->m,input->k * input->d / input->m);
		Cp_index = NE_productQuantA(res,input->n,input->d,input->m,input->k,Cp,input->eps,input->tmin,input->tmax,input->nr);
	

	
		
		

		printf("Creazione della Inverted List\n");
		
		IL= (int**) get_block (sizeof(int*),input->kc);
		bucket=alloc_vector(input->kc);
		

		len_IL=alloc_vector(input->kc);

		for(i=0; i < input->kc; i++){
			bucket[i]=0;
			
			for(j=0;j< input->n;j++){
				if(Cc_index[j]==i)
					bucket[i]++;
					
			}
			IL[i]=alloc_vector(bucket[i]*nodo);
			len_IL[i]=bucket[i];
			bucket[i]=0;
			

		}

		printf("Popolazione dell'Inverted List\n");
		int ind;
		int sjump,sbucket;
		int* L_i;
		for(i=0;i<input->n;i++){
			
			ind = Cc_index[i];
			L_i = IL[ind];
			sbucket = bucket[ind]*nodo;
			
			L_i[sbucket]=i;
			
			

			for(int j = 1; j < input->m+1; j++)
			{

				L_i[sbucket + j] = Cp_index[(j-1)*input->n+i];
				
				
			}
			
			bucket[ind]++;


			
		}

		dealloc_matrix(res);
		printf("Fine index\n");


	}
	else if(input->exaustive == 0 && nmod4==false){


	

		
		Cc= alloc_matrix(input->kc,input->d);
		Cc_index = alloc_vector(input->n);
		NE_k_means_colU(input->ds,input->n,input->d,input->kc,input->eps,Cc_index,Cc,input->tmin,input->tmax,input->nr);
		
		
		printf("Calcolo dei residui\n");
		
		MATRIX res= residuals_col(input->ds,Cc,Cc_index,input->n,input->d);

		
		printf("Quantizzazione dei residui\n");
	

		Cp = alloc_matrix(input->m,input->k * input->d / input->m);
		Cp_index = NE_productQuantU(res,input->n,input->d,input->m,input->k,Cp,input->eps,input->tmin,input->tmax,input->nr);
	

	
		

		printf("Creazione della Inverted List\n");
		
		IL= (int**) get_block (sizeof(int*),input->kc);
		bucket=alloc_vector(input->kc);

		len_IL=alloc_vector(input->kc);

		
		for(i=0; i < input->kc; i++){
			bucket[i]=0;
			for(j=0;j< input->n;j++){
				if(Cc_index[j]==i)
					bucket[i]++;
					
			}
			IL[i]=alloc_vector(bucket[i]*nodo);
			len_IL[i]=bucket[i];
			bucket[i]=0;
			

		}

		printf("Popolazione dell'Inverted List\n");
		int ind;
		int sjump,sbucket;
		int* L_i;
		for(i=0;i<input->n;i++){
			
			ind = Cc_index[i];
			L_i = IL[ind];
			sbucket = bucket[ind]*nodo;
			
			L_i[sbucket]=i;
			
			

			for(int j = 1; j < input->m+1; j++)
			{

				L_i[sbucket + j] = Cp_index[(j-1)*input->n+i];
				
			
				
			}
			
			
			bucket[ind]++;


			
		}


		dealloc_matrix(res);
		printf("Fine index\n");


	}
	else if(input->exaustive == 1 && nmod4==true){
		
		centroids = alloc_matrix(input->m,input->k * input->d / input->m);
		pq = productQuantA(input->ds, input->n, input->d, input->m, input->k, centroids, input->eps, input->tmin, input->tmax);
		
	}
	
	else if(input->exaustive == 1 && nmod4==false){
		centroids = alloc_matrix(input->m,input->k * input->d / input->m);
		pq = productQuantU(input->ds, input->n, input->d, input->m, input->k, centroids, input->eps, input->tmin, input->tmax);
		

	}
	
}


void pqnn_search(params* input) {
	if(input->exaustive==0 && input->symmetric==1 && submod4 == true && dmod4 == true){
		stored_distance=pre_sdcA(Cp,input->sub,input->m,input->k);
		int i,i_w,ind,result,sjump,sbucket;
		int * k_nn;
		float * result_dist;
		float tmp,nn_dis;
		int C_i, z, t;
		int* L_i;
		c_x=alloc_vector(input->m);
		for(i=0;i< input->nq;i++){
			
			int* label_w = w_near_centroidsA(&x_query[i*input->d],Cc,input->kc,input->d,input->w);
		
			float* res_x= residuals_x(&x_query[i*input->d],Cc,label_w,input->w,input->d);
		
			
			k_nn=alloc_vector(input->knn);
			result_dist=alloc_matrix(input->knn,1);
			nn_dis = FLT_MAX;//DBL_MAX;
			
		
			for(i_w = 0 ; i_w < input->w ; i_w++){
				int k_2 = input->k*input->k;
					for(int j=0;j<input->m;j++){
					uj_x = Uj_x( &res_x[i_w*input->d], j, input->m,1,input->d);
							c_x[j] = centXA(&Cp[j*input->sub*input->k], uj_x, input->k, input->sub);
					dealloc_matrix(uj_x);
					
				}	
				
				C_i= label_w[i_w];
				L_i = IL[C_i];

				
				sbucket=bucket[C_i];
			
				for( ind = 0; ind < sbucket; ind++)
				{	
					

					tmp=0;
					for(z=0; z < input->m; z++){
					
						tmp+= stored_distance[z*k_2+c_x[z]*input->k+L_i[ind*nodo+1+z]];
					}

					if(tmp < nn_dis){							
						max_heap(k_nn,result_dist,L_i[ind*nodo],tmp,nn_dis,input->knn,false, &nn_dis);
						
					}
				}
			}
			
			if(c_max_heap == input->knn)
				for(int k=0; k<input->knn; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
			else
			{
				for(int k=0; k<c_max_heap; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
				for(int k=c_max_heap; k<input->knn; k++){
					input->ANN[i*input->knn+k]=-1;
				}
			}

			
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			dealloc_vector(label_w);
			dealloc_matrix(res_x);
			c_max_heap=0;
			pre_max_heap=0;
			
		}
		dealloc_matrix(stored_distance);
		dealloc_vector(L_i);
		dealloc_vector(c_x);
	}
	else if(input->exaustive==0 && input->symmetric==1 && submod4 == true && dmod4 == false){
		stored_distance=pre_sdcA(Cp,input->sub,input->m,input->k);
	
		int i,i_w,ind,result,sjump,sbucket;
		int * k_nn;
		float * result_dist;
		float tmp,nn_dis;
		int C_i, z, t;
		int* L_i;
		c_x=alloc_vector(input->m);
		for(i=0;i< input->nq;i++){
			
			int* label_w = w_near_centroidsU(&x_query[i*input->d],Cc,input->kc,input->d,input->w);
		
			float* res_x= residuals_x(&x_query[i*input->d],Cc,label_w,input->w,input->d);
		
			
			k_nn=alloc_vector(input->knn);
			result_dist=alloc_matrix(input->knn,1);
			nn_dis = FLT_MAX;//DBL_MAX;

		
			for(i_w = 0 ; i_w < input->w ; i_w++){
				int k_2 = input->k*input->k;
					for(int j=0;j<input->m;j++){
					uj_x = Uj_x( &res_x[i_w*input->d], j, input->m,1,input->d);
							printf("-----------------------\n");
					c_x[j] = centXA(&Cp[j*input->sub*input->k], uj_x, input->k, input->sub);
					
					dealloc_matrix(uj_x);
				}	
			
				C_i= label_w[i_w];
				L_i = IL[C_i];

				
				sbucket=bucket[C_i];
			
				for( ind = 0; ind < sbucket; ind++)
				{	
					

					tmp=0;
					for(z=0; z < input->m; z++){
					
						tmp+= stored_distance[z*k_2+c_x[z]*input->k+L_i[ind*nodo+1+z]];
					}

					if(tmp < nn_dis){							
						max_heap(k_nn,result_dist,L_i[ind*nodo],tmp,nn_dis,input->knn,false,&nn_dis);
						
					}
				}
			}
			
			if(c_max_heap == input->knn)
				for(int k=0; k<input->knn; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
			else
			{
				for(int k=0; k<c_max_heap; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
				for(int k=c_max_heap; k<input->knn; k++){
					input->ANN[i*input->knn+k]=-1;
				}
			}

			
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			dealloc_vector(label_w);
			dealloc_matrix(res_x);
			c_max_heap=0;
			pre_max_heap=0;
			
		}
		dealloc_matrix(stored_distance);
		dealloc_vector(L_i);
		dealloc_vector(c_x);
	}
	else if(input->exaustive==0 && input->symmetric==1 && submod4 == false && dmod4 == true){
		stored_distance=pre_sdcU(Cp,input->sub,input->m,input->k);
	
		int i,i_w,ind,result,sjump,sbucket;
		int * k_nn;
		float * result_dist;
		float tmp,nn_dis;
		int C_i, z, t;
		int* L_i;
		c_x=alloc_vector(input->m);
		for(i=0;i< input->nq;i++){
			
			int* label_w = w_near_centroidsA(&x_query[i*input->d],Cc,input->kc,input->d,input->w);
		
			float* res_x= residuals_x(&x_query[i*input->d],Cc,label_w,input->w,input->d);
		
			
			k_nn=alloc_vector(input->knn);
			result_dist=alloc_matrix(input->knn,1);
			nn_dis = FLT_MAX;//DBL_MAX;
		
			for(i_w = 0 ; i_w < input->w ; i_w++){
				int k_2 = input->k*input->k;
					for(int j=0;j<input->m;j++){
					uj_x = Uj_x( &res_x[i_w*input->d], j, input->m,1,input->d);
							c_x[j] = centXU(&Cp[j*input->sub*input->k], uj_x, input->k, input->sub);
					
					dealloc_matrix(uj_x);
				}  
				
				C_i= label_w[i_w];
				L_i = IL[C_i];

				
				sbucket=bucket[C_i];
			
				for( ind = 0; ind < sbucket; ind++)
				{  
				

				tmp=0;
				for(z=0; z < input->m; z++){
				
					tmp+= stored_distance[z*k_2+c_x[z]*input->k+L_i[ind*nodo+1+z]];
				}
				if(tmp < nn_dis){              
					max_heap(k_nn,result_dist,L_i[ind*nodo],tmp,nn_dis,input->knn,false, &nn_dis);
					
				}
				}
			}
			
			if(c_max_heap == input->knn)
				for(int k=0; k<input->knn; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
			else
			{
				for(int k=0; k<c_max_heap; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
				for(int k=c_max_heap; k<input->knn; k++){
					input->ANN[i*input->knn+k]=-1;
				}
			}

			
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			dealloc_vector(label_w);
			dealloc_matrix(res_x);
			c_max_heap=0;
			pre_max_heap=0;
			
		}
		dealloc_matrix(stored_distance);
		dealloc_vector(L_i);
		dealloc_vector(c_x);
	}

	else if(input->exaustive==0 && input->symmetric==1 && submod4 == false && dmod4 == false){
		stored_distance=pre_sdcU(Cp,input->sub,input->m,input->k);
	
		int i,i_w,ind,result,sjump,sbucket;
		int * k_nn;
		float * result_dist;
		float tmp,nn_dis;
		int C_i, z, t;
		int* L_i;
		c_x=alloc_vector(input->m);
		for(i=0;i< input->nq;i++){
			
			int* label_w = w_near_centroidsU(&x_query[i*input->d],Cc,input->kc,input->d,input->w);
			
			float* res_x= residuals_x(&x_query[i*input->d],Cc,label_w,input->w,input->d);
		
			
			k_nn=alloc_vector(input->knn);
			result_dist=alloc_matrix(input->knn,1);
			nn_dis = FLT_MAX;//DBL_MAX;

		
			for(i_w = 0 ; i_w < input->w ; i_w++){
				int k_2 = input->k*input->k;
					for(int j=0;j<input->m;j++){
					uj_x = Uj_x( &res_x[i_w*input->d], j, input->m,1,input->d);
							c_x[j] = centXU(&Cp[j*input->sub*input->k], uj_x, input->k, input->sub);
					
					dealloc_matrix(uj_x);
				}	
				
				C_i= label_w[i_w];
				L_i = IL[C_i];

				
				sbucket=bucket[C_i];
			
				for( ind = 0; ind < sbucket; ind++)
				{	
					

					tmp=0;
					for(z=0; z < input->m; z++){
					
						tmp+= stored_distance[z*k_2+c_x[z]*input->k+L_i[ind*nodo+1+z]];
					}

					if(tmp < nn_dis){							
						max_heap(k_nn,result_dist,L_i[ind*nodo],tmp,nn_dis,input->knn,false, &nn_dis);
						
					}
				}
			}
			
			if(c_max_heap == input->knn)
				for(int k=0; k<input->knn; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
			else
			{
				for(int k=0; k<c_max_heap; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
				for(int k=c_max_heap; k<input->knn; k++){
					input->ANN[i*input->knn+k]=-1;
				}
			}

			
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			dealloc_vector(label_w);
			dealloc_matrix(res_x);
			c_max_heap=0;
			pre_max_heap=0;
			
		}
		dealloc_matrix(stored_distance);
		dealloc_vector(L_i);
		dealloc_vector(c_x);
	}
	else if(input->exaustive==0 && input->symmetric==0 && submod4 == true && dmod4 == true){
		
		int i,i_w,ind,result,sjump,sbucket;
		int * k_nn;
		float * result_dist;
		float tmp,nn_dis;
		int C_i, z, t;
		int* L_i;
		for(i=0;i< input->nq;i++){
			int* label_w = w_near_centroidsA(&x_query[i*input->d],Cc,input->kc,input->d,input->w);
		
			float* res_x= residuals_x(&x_query[i*input->d],Cc,label_w,input->w,input->d);
		
			
			k_nn=alloc_vector(input->knn);
			result_dist=alloc_matrix(input->knn,1);
			nn_dis = FLT_MAX;//DBL_MAX;
		
			for(i_w = 0 ; i_w < input->w ; i_w++){		
				
				stored_distance=pre_adcA(&res_x[i_w*input->d],Cp,input->d,input->m,input->k,input->sub);
				
				
				C_i= label_w[i_w];
				L_i = IL[C_i];
				sbucket=bucket[C_i];

				
				for( ind = 0; ind < sbucket; ind++)
				{
					tmp = 0;
					for(int j=0; j < input->m; j++){
						tmp+=stored_distance[j*input->k+L_i[ind*nodo+1+j]];
					}
					
		
					if(tmp < nn_dis){
						max_heap(k_nn,result_dist, L_i[ind*nodo],tmp,nn_dis,input->knn,false,&nn_dis);
						
					}
					
				}		
		
				dealloc_matrix(stored_distance);
				
			}

			
			if(c_max_heap == input->knn)
				for(int k=0; k<input->knn; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
			else
			{
				for(int k=0; k<c_max_heap; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
				for(int k=c_max_heap; k<input->knn; k++){
					input->ANN[i*input->knn+k]=-1;
				}
			}

			
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			dealloc_matrix(res_x);
			dealloc_vector(label_w);
			

			c_max_heap=0;
			pre_max_heap=0;
		}
	}
	else if(input->exaustive==0 && input->symmetric==0 && submod4 == true && dmod4 == false){
		int i,i_w,ind,result,sjump,sbucket;
		int * k_nn;
		float * result_dist;
		float tmp,nn_dis;
		int C_i, z, t;
		int* L_i;
		for(i=0;i< input->nq;i++){
			int* label_w = w_near_centroidsU(&x_query[i*input->d],Cc,input->kc,input->d,input->w);
		
			float* res_x= residuals_x(&x_query[i*input->d],Cc,label_w,input->w,input->d);
		
			
			k_nn=alloc_vector(input->knn);
			result_dist=alloc_matrix(input->knn,1);
			nn_dis = FLT_MAX;
		
			for(i_w = 0 ; i_w < input->w ; i_w++){		
				
				stored_distance=pre_adcA(&res_x[i_w*input->d],Cp,input->d,input->m,input->k,input->sub);

				
				C_i= label_w[i_w];
				L_i = IL[C_i];
				sbucket=bucket[C_i];


				for( ind = 0; ind < sbucket; ind++)
				{
					tmp = 0;
					for(int j=0; j < input->m; j++){
						tmp+=stored_distance[j*input->k+L_i[ind*nodo+1+j]];
					}
					
		
					if(tmp < nn_dis){
						max_heap(k_nn,result_dist, L_i[ind*nodo],tmp,nn_dis,input->knn,false, &nn_dis);
						
					}
					
				}	
				dealloc_matrix(stored_distance);	
				
			}
			
			if(c_max_heap == input->knn)
				for(int k=0; k<input->knn; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
			else
			{
				for(int k=0; k<c_max_heap; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
				for(int k=c_max_heap; k<input->knn; k++){
					input->ANN[i*input->knn+k]=-1;
				}
			}

			
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			dealloc_matrix(res_x);
			dealloc_vector(label_w);
			c_max_heap=0;
			pre_max_heap=0;
		}
	}
	else if(input->exaustive==0 && input->symmetric==0 && submod4 == false && dmod4 == true){
	
		int i,i_w,ind,result,sjump,sbucket;
		int * k_nn;
		float * result_dist;
		float tmp,nn_dis;
		int C_i, z, t;
		int* L_i;
		for(i=0;i< input->nq;i++){
			int* label_w = w_near_centroidsA(&x_query[i*input->d],Cc,input->kc,input->d,input->w);
		
			float* res_x= residuals_x(&x_query[i*input->d],Cc,label_w,input->w,input->d);
		
			
			k_nn=alloc_vector(input->knn);
			result_dist=alloc_matrix(input->knn,1);
			nn_dis = FLT_MAX;//DBL_MAX;
		
			for(i_w = 0 ; i_w < input->w ; i_w++){		
				
				stored_distance=pre_adcU(&res_x[i_w*input->d],Cp,input->d,input->m,input->k,input->sub);

				
				C_i= label_w[i_w];
				L_i = IL[C_i];
				sbucket=bucket[C_i];


				for( ind = 0; ind < sbucket; ind++)
				{
					tmp = 0;
					for(int j=0; j < input->m; j++){
						tmp+=stored_distance[j*input->k+L_i[ind*nodo+1+j]];
					}
					
		
					if(tmp < nn_dis){
						max_heap(k_nn,result_dist, L_i[ind*nodo],tmp,nn_dis,input->knn,false, &nn_dis);
						
					}
					
				}	
				dealloc_matrix(stored_distance);	
				
			}
			
			if(c_max_heap == input->knn)
				for(int k=0; k<input->knn; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
			else
			{
				for(int k=0; k<c_max_heap; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
				for(int k=c_max_heap; k<input->knn; k++){
					input->ANN[i*input->knn+k]=-1;
				}
			}

			
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			dealloc_matrix(res_x);
			dealloc_vector(label_w);
			c_max_heap=0;
			pre_max_heap=0;
		}
	}
	else if(input->exaustive==0 && input->symmetric==0 && submod4 == false && dmod4 == false){
		
		int i,i_w,ind,result,sjump,sbucket;
		int * k_nn;
		float * result_dist;
		float tmp,nn_dis;
		int C_i, z, t;
		int* L_i;
		for(i=0;i< input->nq;i++){
			int* label_w = w_near_centroidsU(&x_query[i*input->d],Cc,input->kc,input->d,input->w);
		
			float* res_x= residuals_x(&x_query[i*input->d],Cc,label_w,input->w,input->d);
		
			
			k_nn=alloc_vector(input->knn);
			result_dist=alloc_matrix(input->knn,1);
			nn_dis = FLT_MAX;//DBL_MAX;
			
		
			for(i_w = 0 ; i_w < input->w ; i_w++){		
				
				stored_distance=pre_adcU(&res_x[i_w*input->d],Cp,input->d,input->m,input->k,input->sub);

				
				C_i= label_w[i_w];
				L_i = IL[C_i];
				sbucket=bucket[C_i];


				for( ind = 0; ind < sbucket; ind++)
				{
					tmp = 0;
					for(int j=0; j < input->m; j++){
						tmp+=stored_distance[j*input->k+L_i[ind*nodo+1+j]];
					}
					
		
					if(tmp < nn_dis){
						max_heap(k_nn,result_dist, L_i[ind*nodo],tmp,nn_dis,input->knn,false, &nn_dis);
						
					}
					
				}	
				dealloc_matrix(stored_distance);	
				
			}
			
			if(c_max_heap == input->knn)
				for(int k=0; k<input->knn; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
			else
			{
				for(int k=0; k<c_max_heap; k++){
					input->ANN[i*input->knn+k]=k_nn[k];
				}
				for(int k=c_max_heap; k<input->knn; k++){
					input->ANN[i*input->knn+k]=-1;
				}
			}

			
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			dealloc_matrix(res_x);
			dealloc_vector(label_w);
			c_max_heap=0;
			pre_max_heap=0;
		}
	}


	
	else if(input->exaustive==1 && input->symmetric==1 && submod4 == true){
		stored_distance=pre_sdcA(centroids,input->sub,input->m,input->k);

		int k_2 = input->k * input->k;

		int* k_nn;
		float* result_dist;
		float tmp,nn_dis;
		
		c_x=alloc_vector(input->m);
	
		int x,y,k,i;
		for(x=0; x<input->nq;x++){	
			k_nn = alloc_vector(input->knn);	
			result_dist=alloc_matrix(input->knn,1);
			for(int j=0;j<input->m;j++){
				uj_x = Uj_x( &input->qs[x*input->d], j, input->m,1,input->d);
				c_x[j] = centXA(&centroids[j*input->k*input->sub], uj_x, input->k, input->sub);
				dealloc_matrix(uj_x);
			}	
			
		

			nn_dis = FLT_MAX;//DBL_MAX;
			
			for(y=0; y< input->n; y++){

				tmp=0;
				for(int j=0; j< input->m; j++){
					tmp+= stored_distance[j*k_2+c_x[j]*input->k+pq[j*input->n+y]];
				}
					if(tmp < nn_dis){
					max_heap(k_nn,result_dist,y,tmp,nn_dis,input->knn,false, &nn_dis);
					
				}
			}//for y
			
			
			if(c_max_heap == input->knn)
				for(int i=0; i<input->knn; i++){
					input->ANN[x*input->knn+i]=k_nn[i];
				}
			else
			{
				for(int i=0; i<c_max_heap; i++){
					input->ANN[x*input->knn+i]=k_nn[i];
				}
				for(int i=c_max_heap; i<input->knn; i++){
					input->ANN[x*input->knn+i]=-1;
				}
			}
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			c_max_heap=0;
			pre_max_heap=0;
			
		}

		dealloc_vector(c_x);
		dealloc_matrix(stored_distance);

		
	}
	else if(input->exaustive==1 && input->symmetric==1 && submod4 == false){
		stored_distance=pre_sdcU(centroids,input->sub,input->m,input->k);

		int k_2 = input->k * input->k;
		int* k_nn;
		float* result_dist;
		float tmp,nn_dis;
		c_x=alloc_vector(input->m);
		int x,y,k,i;
		for(x=0; x<input->nq;x++){	
			k_nn = alloc_vector(input->knn);	
			result_dist=alloc_matrix(input->knn,1);
			for(int j=0;j<input->m;j++){
				uj_x = Uj_x( &input->qs[x*input->d], j, input->m,1,input->d);
				c_x[j] = centXU(&centroids[j*input->k*input->sub], uj_x, input->k, input->sub);
				dealloc_matrix(uj_x);
			}	
			
		

			nn_dis = FLT_MAX;//DBL_MAX;
			
			for(y=0; y< input->n; y++){

				tmp=0;
				for(int j=0; j< input->m; j++){
					tmp+= stored_distance[j*k_2+c_x[j]*input->k+pq[j*input->n+y]];
				}
					if(tmp < nn_dis){
					max_heap(k_nn,result_dist,y,tmp,nn_dis,input->knn,false, &nn_dis);
					
				}
			}//for y
			
			
			if(c_max_heap == input->knn)
				for(int i=0; i<input->knn; i++){
					input->ANN[x*input->knn+i]=k_nn[i];
				}
			else
			{
				for(int i=0; i<c_max_heap; i++){
					input->ANN[x*input->knn+i]=k_nn[i];
				}
				for(int i=c_max_heap; i<input->knn; i++){
					input->ANN[x*input->knn+i]=-1;
				}
			}
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			c_max_heap=0;
			pre_max_heap=0;
			
		}
		dealloc_vector(c_x);
		dealloc_matrix(stored_distance);

		
	}
	else if(input->exaustive==1 && input->symmetric==0 && submod4 == true){
		float tmp,nn_dis;
		int x,y,k;
		int* k_nn;
		float* result_dist;
		for(x=0; x<input->nq;x++){
			k_nn = alloc_vector(input->knn);
			result_dist=alloc_matrix(input->knn,1);
			stored_distance=pre_adcA(&input->qs[x*input->d],centroids,input->d,input->m,input->k,input->sub);
		

			nn_dis = FLT_MAX;
			for(y=0; y< input->n; y++){
				tmp = 0;
					for(int j=0; j < input->m; j++){
					tmp+=stored_distance[j*input->k+pq[j*input->n+y]];
				}
			
					if(tmp < nn_dis){
					max_heap(k_nn,result_dist,y,tmp,nn_dis,input->knn,false, &nn_dis);
				}
				
			}
			dealloc_matrix(stored_distance);
		


			if(c_max_heap == input->knn)
				for(int i=0; i<input->knn; i++){
					input->ANN[x*input->knn+i]=k_nn[i];
				}
			else
			{
				for(int i=0; i<c_max_heap; i++){
					input->ANN[x*input->knn+i]=k_nn[i];
				}
				for(int i=c_max_heap; i<input->knn; i++){
					input->ANN[x*input->knn+i]=-1;
				}
			}
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			c_max_heap=0;
			pre_max_heap = 0;	
			
		}
		
	}
	else if(input->exaustive==1 && input->symmetric==0 && submod4 == false){
		float tmp,nn_dis;
		int x,y,k;
		int* k_nn;
		float* result_dist;
		for(x=0; x<input->nq;x++){
			
			k_nn = alloc_vector(input->knn);
			result_dist=alloc_matrix(input->knn,1);

			stored_distance=pre_adcU(&input->qs[x*input->d],centroids,input->d,input->m,input->k,input->sub);
		


			nn_dis = FLT_MAX;
			for(y=0; y< input->n; y++){
				tmp = 0;
					for(int j=0; j < input->m; j++){
					tmp+=stored_distance[j*input->k+pq[j*input->n+y]];
				}
			
				
					if(tmp < nn_dis){
					max_heap(k_nn,result_dist,y,tmp,nn_dis,input->knn,false, &nn_dis);
				}
				
			}
			dealloc_matrix(stored_distance);
		


		
			if(c_max_heap == input->knn)
				for(int i=0; i<input->knn; i++){
					input->ANN[x*input->knn+i]=k_nn[i];
				}
			else
			{
				for(int i=0; i<c_max_heap; i++){
					input->ANN[x*input->knn+i]=k_nn[i];
				}
				for(int i=c_max_heap; i<input->knn; i++){
					input->ANN[x*input->knn+i]=-1;
				}
			}
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			c_max_heap=0;
			pre_max_heap = 0;	
			
		}
		
	}


}


int main(int argc, char** argv) {
	char fname[256];
	int i, j;
	
	params* input = malloc(sizeof(params));

	
	input->filename = NULL;
	input->exaustive = 1;
	input->symmetric = 1;
	input->knn = 4;
	input->m = 5;
	input->k = 256;
	input->kc = 256;
	input->w=8;
	input->eps = 0.01;
	input->tmin = 10;
	input->tmax = 100;
	input->silent = 0;
	input->display = 0;
	

	int par = 1;
	while (par < argc) {
		if (par == 1) {
			input->filename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-knn") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing knn value!\n");
				exit(1);
			}
			input->knn = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-m") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing m value!\n");
				exit(1);
			}
			input->m = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-k") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing k value!\n");
				exit(1);
			}
			input->k = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-kc") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing kc value!\n");
				exit(1);
			}
			input->kc = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-w") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing w value!\n");
				exit(1);
			}
			input->w = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-nr") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing nr value!\n");
				exit(1);
			}
			input->nr = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-eps") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing eps value!\n");
				exit(1);
			}
			input->eps = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-tmin") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing tmin value!\n");
				exit(1);
			}
			input->tmin = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-tmax") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing tmax value!\n");
				exit(1);
			}
			input->tmax = atoi(argv[par]);
			par++;
 		} else if (strcmp(argv[par],"-exaustive") == 0) {
 			input->exaustive = 1;
 			par++;
 		} else if (strcmp(argv[par],"-noexaustive") == 0) {
 			input->exaustive = 0;
 			par++;
 		} else if (strcmp(argv[par],"-sdc") == 0) {
 			input->symmetric = 1;
 			par++;
 		} else if (strcmp(argv[par],"-adc") == 0) {
 			input->symmetric = 0;
 			par++;
		} else
			par++;
	}
	
	
	if (!input->silent) {
		printf("Usage: %s <data_name> [-d][-s][-exaustive|-noexaustive][-sdc|-adc][...]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\t-d : display ANNs\n");
		printf("\t-s : silent\n");
		printf("\t-m: PQ groups\n");
		printf("\t-k: PQ centroids\n");
		printf("\t-kc: coarse VQ centroids\n");
		printf("\t-w: coarse VQ centroids to be selected\n");
		printf("\t-nr: residual sample size\n");
		printf("\t-eps: k-means termination threshold\n");
		printf("\t-tmin: min k-means iterations\n");
		printf("\t-tmax: max k-means iterations\n");
		printf("\n");
	}

	if (input->filename == NULL || strlen(input->filename) == 0) {
		printf("Missing input file name!\n");
		exit(1);
	}
	
	sprintf(fname, "%s.ds", input->filename);
	input->ds = load_data_col(fname, &input->n, &input->d);
	input->sub=input->d/input->m;
	

	if(input->nr <= 0 || input->nr > input->n)
		input->nr = input->n/20;


	sprintf(fname, "%s.qs", input->filename);
	input->qs = load_data_row(fname, &input->nq, &input->d);
	

	

	int nmodul= input->n % 8;
	int dmodul= input->d % 8;
	int submodul = (input->d/input->m) % 8;

	if(nmodul == 0)
		nmod4=true;
	if(dmodul == 0)
		dmod4=true;
	if(submodul == 0)
		submod4=true;

	if (input->d % input->m != 0) {
		printf(" Dimensions are not divisible by m \n");
		exit(1);
	}
	if (input->k > input->n || input->kc > input->n) {
		printf(" A lot of centroids compared on dataset \n");
		exit(1);
	}

	if (input->exaustive==0 && input->k > input->nr || input->kc > input->nr) {
		printf(" A lot of centroids compared on dataset \n");
		exit(1);
	}
	
	
	if (!input->silent) {
		printf("Input file name: '%s'\n", input->filename);
		printf("Data set size [n]: %d\n", input->n);
		printf("Number of dimensions [d]: %d\n", input->d);
		printf("Query set size [nq]: %d\n", input->nq);
		printf("Number of ANN [knn]: %d\n", input->knn);
		printf("PQ groups [m]: %d\n", input->m);
		printf("PQ centroids [k]: %d\n", input->k);
		if (!input->exaustive) {
			printf("Coarse VQ centroids [kc]: %d\n", input->kc);
			printf("Coarse VQ centroids to be selected [w]: %d\n", input->w);
			printf("Number of residuals used to determine PQ centroids [nr]: %d\n", input->nr);
		}
		printf("K-means parameters: eps = %.4f, tmin = %d, tmax = %d\n", input->eps, input->tmin, input->tmax);
	}
	
	
	clock_t t = clock();
	pqnn_index(input);
	t = clock() - t;

	
	if (!input->silent)
		printf("\nIndexing time = %.3f secs\n", ((float)t)/CLOCKS_PER_SEC);
	else
		printf("%.3f\n", ((float)t)/CLOCKS_PER_SEC);

	input->ANN = calloc(input->nq*input->knn,sizeof(int));

	clock_t t_1 = clock();
	pqnn_search(input);
	t_1 = clock() - t_1;
	
	if (!input->silent)
		printf("\nSearching time = %.3f secs\n", ((float)t_1)/CLOCKS_PER_SEC);
	else
		printf("%.3f\n", ((float)t_1)/CLOCKS_PER_SEC);
	
	
	
 	if (input->ANN != NULL)
 	{
 		if (!input->silent && input->display) {
 			printf("\nANN:\n");
 			for (i = 0; i < input->nq; i++) {
				printf("\nquery #%d:", i);
				for (j = 0; j < input->knn; j++)
					printf("  %d  ", input->ANN[i*input->knn+j]);
				printf("\n");
 			}
 		}
 		save_ANN(input->filename, input->ANN, input->nq, input->knn);
 	}
	
	if (!input->silent){
		printf("\nDone.\n");
	}

	printf("\nIndexing time = %.3f secs\n", ((float)t)/CLOCKS_PER_SEC);
	printf("\nSearching time = %.3f secs\n", ((float)t_1)/CLOCKS_PER_SEC);

	return 0;
}