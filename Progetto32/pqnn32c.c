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
//#include "kmeans2.c"

#define 	MATRIX		float*
#define 	x_query 	input->qs
#define     nodo		(input->m+1)

typedef struct {
/*
	input->knn = 1;
	input->m = 8;
	input->k = 256;
	input->kc = 8192;
	input->w = 16;
	input->eps = 0.01;
	input->tmin = 10;
	input->tmax = 100;
	input->silent = 0;
	input->display = 0;
*/

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
	//
	// Inserire qui i campi necessari a memorizzare i Quantizzatori
	//
	// ...
	// ...
	// ...
	//
} params;

//variabili utili per l'algoritmo esaustivo
float * centroids;
int * pq; 
float* uj_x;
int* c_x;
int c_max_heap=0;
float pre_max_heap=0;


//variabili utili per l'algoritmo non esaustivo
MATRIX Cc;
int* Cc_index;
float* Cp;
int * Cp_index;
float* stored_distance;
int* len_IL;
int ** IL;
int* bucket;
//int* jump;


//test
int accesso=0; //utilizzato per SDC
int accesso_1=0; //utilizzato per ADC
int accesso_2=0; //utilizzato per gli accessi nel max_heap
int accesso_3=0;




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
	return _mm_malloc(elements*size,16); 
}


void free_block(void* p) { 
	_mm_free(p);
}


MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(float),rows*cols);
}

//riguardare bene il fatto dell'allineamento (di default è 16)
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

MATRIX load_data(char* filename, int *n, int *d) {	
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
void printCentroids(MATRIX C,int* index, int n,int d,int k){/////////////////////////////////////forse va cambiato///////////////////
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
			printf("ds[%d][%d]=   %f   \n",i,j,ds[i*d+j]);
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
}

//print per testare il metodo Uj
void printEq(MATRIX m1, MATRIX m2, int m1_n,int m1_d,int m2_n,int m2_d){
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

/**metodo per estrapolare in maniera semi-casuale nr elementi da
 * un dataset */

MATRIX extrac(MATRIX ds,int n,int d,int nr){
	MATRIX result = alloc_matrix(nr,d);
	int i,j,h;
	for (h = i = 0; i < nr; h += n/nr , i++) {
		for (j = 0; j < d;j++){
			result[i*d+j] = ds[h*d+j];
		}
	}
	return result;
	
}

//metodo che prende un sottogruppo (sub dimensionale) del data set
// j serve per prendere il j-esimo gruppetto di dimensioni j=2 equivale ad U2
MATRIX Uj(MATRIX ds, int j,int m,int n,int d){
	int i,k,c;
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

float dist(float * x,float * y, int d){
	float distance = 0;
	for (int i=0; i<d;i++){
	    distance += pow(x[i] - y[i], 2);
	}
	return distance;
}

int centX(float * centroids, float * x, int k, int d){	
	float dis = dist(x, centroids, d);
	int park = 0;
	for(int i=1; i<k; i++){
		float tmp = dist(x, &centroids[i*d], d);	
 		if( tmp < dis){
			dis = tmp;
			park = i;
		}
	}
	return park;
}
/*
int centX(float * centroids, float * x, int k, int d){	
	float dis = 0;
	for (int i=0; i<d;i++){
	    dis += pow(x[i] - centroids[i], 2);
	}
	int park = 0;
	for(int i=1; i<k; i++){
		float tmp = 0;
		for (int j=0; j<d;j++){
			tmp += pow(x[j] - centroids[i*d+j], 2);
		}	
 		if( tmp < dis){
			dis = tmp;
			park = i;
		}
	}
	return park;
}*/

/**MATRIX randCentroid(MATRIX ds,int n,int d,int k){
	int i,j;
	int max=0;//rand() % (n-k);
	printf("---------%d-----------",max);
	MATRIX initialCentroid = alloc_matrix(k,d);
	for(i = max; i < max+k; i++){
		for(j = 0; j < d; j++){
			initialCentroid[i*d+j]=ds[i*d+j];
		}
	}
	return initialCentroid;
}**/


int mapping(int i,int j,int n){
	if(i == j ) 
		return 0;
	else if ( i > j ){
		return (n*(n-1)/2) - (n-j)*((n-j)-1)/2 + i - j - 1;
		/*int tmp=i;
		i=j;
		j=tmp;
		//free(&tmp); scoprire se ha senso farlo o no sulle variabili*/
	}else
	{
		return (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
	}
	

	//si può ottimizzare salvando una variabile dim
	//che cmq serve anche per allocare le matrici
	
}

void interClusterCalc(float* MCD, float* centroids,float* stored_distance,int n,int d,int k){
	int i,j,ind;
	float tmp;
	for(i = 0; i < k; i++){
		tmp=0;
		for(j = i+1; j < k; j++){
			tmp = dist(&centroids[i*d],&centroids[j*d],d);
			ind=mapping(i,j,k);
			stored_distance[ind]=tmp;
			if(MCD[i]==-1 || MCD[i] > 0.5*tmp){
				MCD[i]==0.5*tmp;
			}
			if(MCD[j]==-1 || MCD[j] > 0.5*tmp){
				MCD[j]==0.5*tmp;
			}
		}
	}
}





//kmeans modificato in modo da prendere due "MATRIX" e usando l'alloc del prof con l'allineamento.
void k_means(MATRIX data, int n, int d, int k, float t, int* labels, MATRIX centroids,int t_min,int t_max) {
	
	/* output cluster label for each data point */
	//int * labels = alloc_vector(n);

	//queste variabili sono da liberare alla fine del metodo!!!
	float min_distance;
	float distance;
	float offset;
	int iter=0;
	int h, i, j; /* loop counters, of course */
	
	/* size of each cluster */
	int* counts = alloc_vector(k);
	float old_error, error = FLT_MAX;//DBL_MAX; /* sum of squared euclidean distance */
	
	MATRIX c = centroids;
	
	/* temp centroids */
	MATRIX c1 = alloc_matrix(k,d);
	
	/****
	 ** initialization */
	printf("initilization\n");
	for (h = i = 0; i < k; h += n / k, i++) {
		/* pick k points as initial centroids */
		//printf("----h=%d-----i=%d-------\n",h,i);
		for (j = 0; j < d;j++){
			c[i*d+j] = data[h*d+j];
			//printf("c[%d][%d] = data[%d][%d] ------>  %f  ||||  %f \n",i,j,h,j,c[i*d+j],data[h*d+j]);
		}
	}
	
	printf("main loop\n");

	do {
		iter++;
		/* save error from last step */
		old_error = error, error = 0;
		/* clear old counts and temp centroids */
		for (i = 0; i < k;i++) {
			counts[i] = 0;
			for (j = 0; j < d; j++){
				c1[i*d+j] = 0;
			}
		}
		printf("identify the closest cluster in %d iteration\n",iter);
		for (h = 0; h < n; h++) {
			/* identify the closest cluster */
			min_distance = FLT_MAX;//DBL_MAX;
			for (i = 0; i < k; i++) {
				distance = 0;
				for (j = 0; j < d ; j++){
					distance += pow(data[h*d+j] - c[i*d+j], 2);
				}
				if (distance < min_distance) {
					labels[h] = i;
					min_distance = distance;
				}
			}
			//printf("update size and temp centroid of the destination cluster for %d point\n",h);
			/* update size and temp centroid of the destination cluster */
			for (j = 0; j < d; j++){
				c1[labels[h]*d+j] += data[h*d+j]; // c'era un +=
			}
			counts[labels[h]]++;
			/* update standard error */
			error += min_distance;
		}
		printf("Update all centroids\n");
		for (i = 0; i < k; i++) { /* update all centroids */
			for (j = 0; j < d; j++) {
				if(counts[i]!=0){
					c[i*d+j] = c1[i*d+j] / counts[i];
					//printf("SI ---> c[%d][%d] =  %f \n",i,j,c[i*d+j]);
				}else
				{
					c[i*d+j] =c1[i*d+j];
					//printf("NO ---> c[%d][%d] =  %f \n",i,j,c[i*d+j]);
				}
			}
		}
	}//while (fabs(error-old_error) > t); 
	while (!(t_min <= iter && ((t_max < iter) || fabs(error-old_error) <= t)));

	/****
	 ** housekeeping */

	/*for (i = 0; i < k; i++) {
		if (!centroids) {
			free(c[i]);
		}
		free(c1[i]);
	}*/
	
	printf("housekeeping\n");
	
	dealloc_matrix(c1);

	dealloc_vector(counts);

	//return labels;
}//k_means

void triangle_k_means(MATRIX data,float* stored_distance, int n, int d, int k, float t, int* labels, MATRIX centroids,int t_min,int t_max) {
	

	float tmp;
	float min_distance;
	float distance;
	float offset;
	int iter=0;
	int h, i, j; /* loop counters, of course */
	

	
	/* size of each cluster */
	int* counts = alloc_vector(k);
	float old_error, error = FLT_MAX;//DBL_MAX; /* sum of squared euclidean distance */
	
	MATRIX c = centroids;
	
	/* temp centroids */
	MATRIX c1 = alloc_matrix(k,d);
	
	/****
	 ** initialization */
	printf("initilization\n");
	for (h = i = 0; i < k; h += n / k, i++) {
		/* pick k points as initial centroids */
		//printf("----h=%d-----i=%d-------\n",h,i);
		for (j = 0; j < d;j++){
			c[i*d+j] = data[h*d+j];
			//printf("c[%d][%d] = data[%d][%d] ------>  %f  ||||  %f \n",i,j,h,j,c[i*d+j],data[h*d+j]);
		}
	}

	float * l = calloc(n*k,sizeof(float));
	float* u = calloc(n*k,sizeof(float));
	float* ICD = stored_distance;
	//dovrebbe metterli tutti a 0 e quindi tutti false;
	bool* r=calloc(n,sizeof(bool));
	
	float * MCD = alloc_matrix(k,1);
	for(int i =0;i<k;i++){
		MCD[i]=-1;
	}
	interClusterCalc(MCD,centroids,ICD,n,d,k);

	//puliziatmp
	for (i = 0; i < k;i++) {
		counts[i] = 0;
		for (j = 0; j < d; j++){
			c1[i*d+j] = 0;
		}
	}

	//initilization
	for (int i = 0; i < n; i++)
	{
	/* identify the closest cluster */
		int pos=0;
		distance = dist(&data[i*d],&centroids[pos*d],d);
		min_distance=distance;
		labels[i]=pos;
		int currentCentroid = pos;

		//lowerbound[point.id][closestCentroid]
		l[i*k+pos]=min_distance;

		for (int i_c = 0; i_c < k; i_c++) {
			if (i_c != currentCentroid){
				int ind = mapping(labels[i],pos,k);
				if ( 0.5*ICD[ind] < min_distance)
				{
					distance = dist(&data[i*d],&centroids[i_c*d],d);
					l[i*k+pos]=distance;
					if (min_distance > distance) 
					{
					//labels[i] = pos;
					labels[i]=i_c;
					min_distance = distance;
					}
				}
			}
			//pos++;
		}
		u[i]=min_distance;
		counts[labels[i]]++;

		for (j = 0; j < d; j++){
				c1[labels[i]*d+j] += data[i*d+j]; // c'era un +=
		}
		/* update standard error */
		error += min_distance;
	
	}

	//recalculateCentroids()
	int pos=0;
	printf("Update all centroids\n");
	for (i = 0; i < k; i++) { /* update all centroids */
		for (j = 0; j < d; j++) {
			if(counts[i]!=0){
				c[i*d+j] = c1[i*d+j] / counts[i];
			}else
			{
				c[i*d+j] =c1[i*d+j];
			}
		}
	}

	do {
			iter++;
			/* save error from last step */
			old_error = error, error = 0;
			/* clear old counts and temp centroids */
			for (i = 0; i < k;i++) {
				counts[i] = 0;
				for (j = 0; j < d; j++){
					c1[i*d+j] = 0;
				}
			}

			//assignPoints
			interClusterCalc(MCD,centroids,ICD,n,d,k);
			


			for (h = 0; h < n; h++) {


				/* identify the closest cluster */
				min_distance = FLT_MAX;//DBL_MAX;
				for (i = 0; i < k; i++) {
					distance = 0;
					for (j = 0; j < d ; j++){
						distance += pow(data[h*d+j] - c[i*d+j], 2);
					}
					if (distance < min_distance) {
						labels[h] = i;
						min_distance = distance;
					}
				}
				//printf("update size and temp centroid of the destination cluster for %d point\n",h);
				/* update size and temp centroid of the destination cluster */
				for (j = 0; j < d; j++){
					c1[labels[h]*d+j] += data[h*d+j]; // c'era un +=
				}
				counts[labels[h]]++;
				/* update standard error */
				error += min_distance;
			}
			printf("Update all centroids\n");
			for (i = 0; i < k; i++) { /* update all centroids */
				for (j = 0; j < d; j++) {
					if(counts[i]!=0){
						c[i*d+j] = c1[i*d+j] / counts[i];
						//printf("SI ---> c[%d][%d] =  %f \n",i,j,c[i*d+j]);
					}else
					{
						c[i*d+j] =c1[i*d+j];
						//printf("NO ---> c[%d][%d] =  %f \n",i,j,c[i*d+j]);
					}
				}
			}
		}//while (fabs(error-old_error) > t); 
		while (!(t_min <= iter && ((t_max < iter) || fabs(error-old_error) <= t)));

		
		printf("housekeeping\n");
		
		dealloc_matrix(c1);

		dealloc_vector(counts);

		//return labels;
	}//k_means


MATRIX residuals(MATRIX ds,MATRIX centroids,int* label, int n,int d){
	MATRIX results = alloc_matrix(n,d);
	int i,j;
	for (i=0;i< n;i++){
		for(j=0;j<d;j++){
			results[i*d+j]=ds[i*d+j]-centroids[label[i]*d+j];
		}
	}
	
	return results;
}

/**in questo modo il primo indice è j per individuare il sottogruppo di centroidi e label.
centroids viene inizializzato all'interno quindi va solo passato un puntatore vuoto.
il primo indice indica il gruppetto il secondo invece indica la dimensione.
si è scelto di rimanere coerenti con le altre implementazioni, per ora tutto cioè che viene 
passato come parametro ci si aspetta sia già allocato mentre tutto cioè che sta dentro il metodo
compreso il valore di ritorno si alloca dentro il metodo.*/
int* productQuant(MATRIX ds,int n,int d,int m,int k,float* centroids,float eps,int t_min,int t_max){
	int j;
	int sub=d/m;
	//int** result = (int**) get_block(sizeof(int*),m);
	int* result = alloc_vector(m*n);
	//centroids = (float**) get_block(sizeof(MATRIX),m);
	for( j = 0; j < m; j++)
	{
		//printf("\nCalcolo del %d sotto-gruppo di centroidi\n",j);
		if(m!=1){
			MATRIX tmp = Uj(ds,j,m,n,d);
			//centroids[j]=alloc_matrix(k,sub);//////////////////////////////////////////////////////
			//result[j]=k_means(tmp,n,sub,k,eps,&centroids[j*sub*k],t_min,t_max);
			k_means(tmp,n,sub,k,eps,&result[j*n],&centroids[j*sub*k],t_min,t_max);
			dealloc_matrix(tmp); // da testare
		}else
		{
			k_means(ds,n,sub,k,eps,&result[j*n],&centroids[j*sub*k],t_min,t_max);
		}
		
	}
	return result;
}


/**
 * metodo per mantenere i k elementi più vicini
 * come se fossero in un heap.
 * In questo caso si è scelto di non mantenere ordinata la struttura ma
 * di utilizzare una variabile per mantenere il massimo corrente in modo
 * da utilizzarla come filtro per risparmiare accessi in memoria
 * si cercano di mantenere i nomi utilizzati nell'algortimo w_nearest_centroids
 * se full == true allora la struttura si comporta come previsto se invece la variabile full
 * è false allora il metodo per,finchè non sarà completamente popoalta, carica la struttura
 * con gli elementi e ritorna il massimo 
 **/
float max_heap(int* index,float* result_dist,int y,float tmp,float max,int dim,bool full){
	//accesso_2++;
	float new_max;
	if (full==true || c_max_heap == dim) { // dovrebbe andar bene anche solo con ==
		//printf("IF:\n C=%d\n",c_max_heap);
		bool trovato = false;
		new_max=tmp;
		//bisogna inserire ed aggiornare la struttura
		for(int j = 0; j < dim; j++)
		{
			//ho il dubbio che ci possano essere più punti con la stessa distanza
			//ci può stare un controllo
			if(!trovato && result_dist[j]==max){
				//printf("cambio NN tolgo il punto v[%d] = %f e metto quello %d\n",j,result_dist[j],y);
				result_dist[j]=tmp;
				index[j]=y;
				trovato=true;
				//printVector(result_w,w);
			}else if (result_dist[j] > new_max)
			{
				new_max = result_dist[j];
			}	
		}
		return new_max;
	}
	else
	{
		//printf("ELSE:\n C=%d\n",c_max_heap);
		index[c_max_heap]=y;
		result_dist[c_max_heap]=tmp;
		//printf("index[%d]=%d\n result_dist[%d]=%f\n",c_max_heap,y,c_max_heap,result_dist[c_max_heap]);
		if (tmp > pre_max_heap) {
			pre_max_heap=tmp;
		}
		c_max_heap++;
		if (c_max_heap==dim) {
			return pre_max_heap;
		}else
		{
			return max;
		}
	}
}



/*x=query
n è il numero dei centroidi del quantizzatore coarse
d dimensione dei centroids
w numero di centroidi "vicini" da analizzare*/
/*int * w_near_centroids(MATRIX x,MATRIX centroids,int n,int d,int w){
	int i,j;
	int * result_w=alloc_vector(w);
	float * result_dist=alloc_matrix(w,1);
	float tmp=0;
	float max=0;

	//riempo i primi w posti con i primi w centroidi e le relative distanze
	printf("riempo i primi w posti\n");
	for(i = 0; i < w; i++)
	{	
		tmp=dist(x,&centroids[i],d);
		result_w[i]=i;
		result_dist[i]=tmp;
		//piccola ottimizzazione, al posto di mantenere ordinata la struttura
		//uso un "max" se le distanze che calcolo sono più piccole allora dovrà 
		//entrare nella struttura altrimenti no.
		if (max < tmp) {
			max = tmp;
		}
	}

	int new_i;
	float new_max;
	bool trovato;
	//n qui simboleggia il numero dei centroidi
	
	printf("incomincio a analizzare tutti i centroidi per il calcolo dei w più vicini\n");
	for(i=w;i<n;i++){
		tmp=dist(x,&centroids[i],d);
		//printf("\nil centroide num[%d] con X dista = %f\n",i,tmp);
		//printf("la distanza max della struttura è = %f\n",max);

		if(tmp < max){
			new_max=tmp;
			//bisogna inserire ed aggiornare la struttura
			for(int j = 0; j < w; j++){
				//ho il dubbio che ci possano essere più punti con la stessa distanza
				//ci può stare un controllo
				if(!trovato && result_dist[j]==max){
					//printf("cambio centroide tolgo il centroide v[%d] = %d e metto quello %d\n",j,result_w[j],i);
					result_dist[j]=tmp;
					result_w[j]=i;
					trovato=true;
					//printVector(result_w,w);

				}else if (result_dist[j] > new_max)
				{
					new_max = result_dist[j];
				}
				
			}
			max = new_max;
			trovato=false;
			
		}

	}
	dealloc_matrix(result_dist);
	return result_w;
}*/

int * w_near_centroids(MATRIX x,MATRIX centroids,int n,int d,int w){
	int i,j;
	int * result_w=alloc_vector(w);
	float * result_dist=alloc_matrix(w,1);
	float tmp=0;
	float max=0;

	//riempo i primi w posti con i primi w centroidi e le relative distanze
	//printf("riempo i primi w posti\n");
	for(i = 0; i < w; i++)
	{	
		//tmp=dist(x,&centroids[i],d);
		tmp = 0;
		for (int j=0; j<d;j++){
			tmp += pow(x[j] - centroids[i*d+j], 2);
		}
		result_w[i]=i;
		result_dist[i]=tmp;
		//piccola ottimizzazione, al posto di mantenere ordinata la struttura
		//uso un "max" se le distanze che calcolo sono più piccole allora dovrà 
		//entrare nella struttura altrimenti no.
		if (max < tmp) {
			max = tmp;
		}
	}

	int new_i;
	//	float new_max;
	//	bool trovato;
	//n qui simboleggia il numero dei centroidi
	
	//printf("incomincio a analizzare tutti i centroidi per il calcolo dei w più vicini\n");
	for(i=w;i<n;i++){
		//tmp=dist(x,&centroids[i],d);
		//printf("\nil centroide num[%d] con X dista = %f\n",i,tmp);
		//printf("la distanza max della struttura è = %f\n",max);
		tmp = 0;
		for (int j=0; j<d;j++){
			tmp += pow(x[j] - centroids[i*d+j], 2);
		}
		
		if(tmp < max){
	//			new_max=tmp;
			//bisogna inserire ed aggiornare la struttura
			max=max_heap(result_w,result_dist,i,tmp,max,w,true);
	//			max = new_max;
	//			trovato=false;
			
		}

	}
	dealloc_matrix(result_dist);
	return result_w;


}

//in questo metodo la x si presuppone 1 vettore d dimensionale quindi non una matrice.
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



float sdc(int* c_x,float* stored_distance, int y, int m,int n, int* labels, int k ){
	float dis=0;
	int i,j;
	for(j=0; j< m; j++){
		
		//old_dis += pow(dist(& centroids[j][c_x*d/m],& centroids[j][labels[j][y]*d/m],d/m),2);

		//questo controllo è dovuto al fatto che la funzione mapping ritorna l'indice corretto quando è
		//possibile altrimenti quando i==j torna direttamente 0 e a sto punto evito di accedere alla struttura
		//i=mapping(c_x[j],labels[j][y],k);
		i=mapping(c_x[j],labels[j*n+y],k);
		if (i!=0) {
			//dis+= stored_distance[j][i];
			dis+= stored_distance[j*(k*(k-1)/2)+i];
		}
		//printf("\nold_dis = %f\ndis = %f\n c_x = %d , label[%d][%d] = %d , i=%d\n",old_dis,dis,c_x,j,y,labels[j][y],i);
	}
	//printf("SDC\nold_dis = %f\ndis = %f\n" ,old_dis,dis);
	return dis;
}


//ancora da smaltire
float adc(float* stored_distance, int y, int k, int m, int n, int* labels){
	float dis = 0;
	for(int j=0; j<m; j++){
		//old_dis += pow(dist(uj_x, & centroids[j][labels[j][y]*d/m],d/m),2);
		dis+=stored_distance[j*k+labels[j*n+y]];
	}
	//printf("ADC\nold_dis = %f\ndis = %f\n" ,old_dis,dis);
	return dis;
}


//a differenza di prima qui abbiamo qua tutte le info calcolate anche per quanto riguarda 
//il quantizzato di r(y) e quindi inutile usare la label per ottenere le info ma conviene 
//passarle direttamente in input
float NE_adc(float* stored_distance, int k, int m,int* res){
	float dis = 0;
	for(int j=0; j<m; j++){
		//old_dis += pow(dist(uj_x, & centroids[j][labels[j][y]*d/m],d/m),2);
		dis+=stored_distance[j*k+res[j]];
	}
	//printf("ADC\nold_dis = %f\ndis = %f\n" ,old_dis,dis);
	return dis;
}
float NE_sdc(int* c_x,float* stored_distance,int m, int* res,int k ){
	float dis=0;
	int i,j;
	for(j=0; j< m; j++){
		
		//old_dis += pow(dist(& centroids[j][c_x*d/m],& centroids[j][labels[j][y]*d/m],d/m),2);

		//questo controllo è dovuto al fatto che la funzione mapping ritorna l'indice corretto quando è
		//possibile altrimenti quando i==j torna direttamente 0 e a sto punto evito di accedere alla struttura
		i=mapping(c_x[j],res[j],k);
		//printf(" i=%d  ------- c_x[%d]=%d -------- res[%d]=%d\n",i,j,c_x[j],j,res[j]);
		if (i!=0) {
			dis+= stored_distance[j*(k*(k-1)/2)+i];
		}
		//printf("\nold_dis = %f\ndis = %f\n c_x = %d , label[%d][%d] = %d , i=%d\n",old_dis,dis,c_x,j,y,labels[j][y],i);
	}
	//printf("SDC\nold_dis = %f\ndis = %f\n" ,old_dis,dis);
	return dis;
}

/*
x per noi è 1 solo punto d dimensionale
ottimizzabile come gli altri con una riga in modo
da avanzare nel gruppo di centroidi j
*/

float* pre_adc(MATRIX x, float* centroids,int d,int m, int k ){
	//float* result=(float**)get_block(sizeof(float*),m);
	float* result= alloc_matrix(m,k);
	int sub=d/m;
	int i,j;
	float distance;
	MATRIX uj_x;
	for(j=0; j<m; j++){
		uj_x = Uj( x, j, m, 1, d);
		//result[j]=alloc_matrix(k,1);
		for(i = 0; i < k; i++){
			//result[j*k+i] = dist(uj_x, &centroids[j*k*sub+i*sub],sub);
			distance = 0;
			for (int z=0; z<sub ;z++){
				distance += pow(uj_x[z] - centroids[j*k*sub+i*sub+z], 2);
			}
			result[j*k+i] = distance;
			//printf("\ncalcolo della distanza U_x[%d] e C[%d][%d] = %f\n",j,j,i,result[j][i]);

		}

		//dis += pow(dist(uj_x, & centroids[j][labels[j][y]*d/m],d/m),2);
	}
	return result;
}


/*popolazione con struttura dimezzata delle distanze tra tutti i centroidi di un sottogruppo j
per accedere alla distanza bisogna usare la funzione mapping che ritorna l'indice corretto
trasformando opportunamente gli indici i,j*/
float* pre_sdc(float* centroids,int d,int m, int k ){
	//float** result=(float**)get_block(sizeof(float*),m);
	float* result= alloc_matrix(m,k*(k-1)/2);
	int sub=d/m;
	int i,j,c,j_d;
	float distance;
	for(j=0; j<m; j++){
		//result[j]=alloc_matrix(k*(k-1)/2,1);
		c=0;
		for(i = 0; i < k; i++){
			for(j_d = i+1; j_d < k;j_d++){
				//result[j][c] = dist(&centroids[j*k*sub+i*sub], &centroids[j*k*sub+j_d*sub],sub);
				//result[j*(k*(k-1)/2)+c] = dist(&centroids[j*k*sub+i*sub], &centroids[j*k*sub+j_d*sub],sub);

				
				distance = 0;
				for (int z=0; z<sub;z++){
					distance += pow(centroids[j*k*sub+i*sub+z] - centroids[j*k*sub+j_d*sub+z], 2);
				}

				//funzione NASM
				//distance=rowdistance32(centroids[j*k*sub+i*sub],centroids[j*k*sub+j_d*sub],d);

				result[j*(k*(k-1)/2)+c]=distance;
				//printf("\ncalcolo della distanza C[%d][%d] e C[%d][%d] = %f\n",j,i,j,j_d,result[j][c]);
				c++;
			}
		}

		//dis += pow(dist(uj_x, & centroids[j][labels[j][y]*d/m],d/m),2);
	}
	return result;
}


//extern void pqnn32_index(params* input);
//extern int* pqnn32_search(params* input);
extern void residual_nasm(float* res, float* ds,float* Cc,int* Cc_index, int n, int d);
extern float rowdistance32(float * c1,float* c2,int d);

/*
 *	pqnn_index
 * 	==========
 */
void pqnn_index(params* input) {
	int i,j;

	if(input->exaustive == 0){


		//TEST 
		//era per provare il kmeans1.c
		//MATRIX Cc = randCentroid(input->ds,input->n,input->d,input->kc);
		printDsQs(input->ds,input->qs,input->n,input->d,input->nq);
		//return;


		printf("Quantizzazione y in qc\n");
		//quantizzare y in qc(y) = Ci , prima si crea il "quantizzatore" richiamando k-means
		//qui dovremmo usare un sottoinsieme


		MATRIX sub_set = extrac(input->ds,input->n,input->d,input->nr);
		input->n = input->nr; // per non cambiare tutto dopo

		//Creazione del quantizzatore Coarse e relativi Centroidi

		//printDsQs(input->ds,sub_set,input->n,input->d,input->nr);
		Cc= alloc_matrix(input->kc,input->d);
		Cc_index = alloc_vector(input->n);
		k_means(sub_set,input->n,input->d,input->kc,input->eps,Cc_index,Cc,input->tmin,input->tmax);
		//k_means(input->ds,input->n,input->d,input->kc,input->eps,Cc_index,Cc,input->tmin,input->tmax);
		//printCentroids(Cc,Cc_index,input->n,input->d,input->kc);
		
		printf("Calcolo dei residui\n");
		//calcolo dei redisui r(y) = y - Ci , qui potrei anche usare il sub_set ma secondo
		//me non avrebbe senso quindi per adesso calcoliamo TUTTI i residui per tutto il dataset
		MATRIX res= residuals(sub_set,Cc,Cc_index,input->n,input->d);

		//funzione NASM per il calcolo dei residui (non riesco ad ottimizzarlo in nasm)
		//residuals_nasm(res,input->ds,Cc,Cc_index,input->n,input->d);

		/**in questa variante i residui vengono calcolati su un sottoinsieme di punti del data-set
		 * precedentemenete clusterizzato*/
		
		//MATRIX res= residuals(sub_set,Cc,Cc_index,input->n,input->d);
		//printDsQs(res,NULL,input->n,input->d,0);

		printf("Quantizzazione dei residui\n");
		//quantizzare r(y) con Qp, prima si crea il quantizzatore usando m volte k-means
		//Cp = (float**)get_block(sizeof(float*),input->m);
		Cp = alloc_matrix(input->m,input->k * input->d / input->m);
		Cp_index = productQuant(res,input->n,input->d,input->m,input->k,Cp,input->eps,input->tmin,input->tmax);
		//Cp_index = productQuant(input->ds,input->n,input->d,input->m,input->k,Cp,input->eps,input->tmin,input->tmax);
		/*
		for(int j= 0; j < input->m; j++)
		{
			printCentroids(Cp[j],Cp_index[j],input->n,input->d/input->m,input->k);
		}
		*/
		

		printf("Creazione della Inverted List\n");
		/*aggiungere nella inverted List una tupla corrispondente a qc(y)=Ci
		succesivamente appendere un "oggetto" composto dal l'ID del punto y in analisi
		e l'indice del centroide prodotto sul residuo Qp(res(y))=Cpi */

		//per ora uso la maniera più stupida e creo tutto poi qui sicuramente si può ottimizzare
		//allocazione dinamica
		//IL=(int***) get_block(sizeof(int**),input->kc);
		IL= (int**) get_block (sizeof(int*),input->kc);
		bucket=alloc_vector(input->kc);
		//jump=alloc_vector(input->kc);

		len_IL=alloc_vector(input->kc);

		//sto creando una struttura che ospita sia il numero di centroidi dentro una lista e
		//anche il numero di celle da sorpassare per arrivarci- è importante far avanzare il for di +2
		//jump[0]=0;
		for(i=0; i < input->kc; i++){
			bucket[i]=0;
			//jump[i]=jump[i-1];
			for(j=0;j< input->n;j++){
				if(Cc_index[j]==i)
					//jump[i]++;
					bucket[i]++;
					
			}
			IL[i]=alloc_vector(bucket[i]*nodo);
			len_IL[i]=bucket[i];
			bucket[i]=0;
			

			//printf("jump[%d]=%d\n",i,jump[i]);
			//printf("bucket[%d]=%d\n",i,bucket[i]);
			//IL[i]=(int**)get_block(sizeof(int*),bucket[i]);

			// mi serve questo perchè dopo devo scorrere tutte le celle
			//mentre bucket viene usato solo per riempire la IL e perde le sue info
			//forse possiamo fare a meno di questa struttura len_IL
			//len_IL[i]= bucket[i];

			//fatto per popolare la struttura in maniera coerente alla scansione del dataset
			//bucket[i] = 0;
		}

		printf("Popolazione dell'Inverted List\n");
		int ind;
		int sjump,sbucket;
		int* L_i;
		for(i=0;i<input->n;i++){
			//posso usare bucket[i] per sapere la dim di ogni cosa e magari farmene una copia
			//potrei caricare i valori al contrario facendo -- in modo da essere sicuro che riempo tutto
			//per simulare la tupla vado a creare un vettore dove il primo elemento è l'iD della y e i restanti
			//sono gli indici che compongono il product quantizzazion e quindi m.
			
			//int* nodo = alloc_vector(input->m+1);
			//nodo[0]=i;
			ind = Cc_index[i];
			L_i = IL[ind];
			//sjump = jump[ind]*nodo;
			sbucket = bucket[ind]*nodo;
			
			L_i[sbucket]=i;
			
			//IL[sjump + sbucket]=i;
			

			//printf("\n");
			//printf("Cc_index[%d] = %d\n",i,ind);
			//printf("bucket[%d] = %d\n",ind,bucket[ind]);
			//printf("STAMPA DEL NODO sotto: %d\n",i);
			for(int j = 1; j < input->m+1; j++)
			{

				L_i[sbucket + j] = Cp_index[(j-1)*input->n+i];
				
				//nodo[j]=Cp_index[j-1][i]; // prendo i vari gruppi 
				//forse posso addirittura deallocare Cp_index che avanti non viene usato
				
			}
			
			//bucket[ind]--;
			//IL[ind][bucket[ind]]=nodo;

			//per riuscire ad avanzare del giusto numero di posizioni devi sommarle

			
			bucket[ind]++;


			
		}

		//printVector(bucket,input->kc);


		//bucket dovrebbe essere tutto zero
		//dealloc_vector(bucket);
		printf("Fine index\n");

		if(input->symmetric==1)
		{
			printf("PRE-calcolo delle distanze (SIMMETRICO) tra Cji e Cji\n");
			stored_distance=pre_sdc(Cp,input->d,input->m,input->k);
		}

	}

	if(input->exaustive == 1){
		centroids = alloc_matrix(input->m,input->k * input->d / input->m);
		pq = productQuant(input->ds, input->n, input->d, input->m, input->k, centroids, input->eps, input->tmin, input->tmax);
		//printf("ho calcolato i centroidi (productQuant)\n");

		
		if(input->symmetric==1){
			printf("PRE-Calcolo le distanze (SIMMETRICO) tra Cji e Cji\n");
			stored_distance=pre_sdc(centroids,input->d,input->m,input->k);
		}


	}
    // -------------------------------------------------
    // Codificare qui l'algoritmo di indicizzazione
    // -------------------------------------------------
    
    //pqnn32_index(input); // Chiamata funzione assembly

    // -------------------------------------------------

}


/*
 *	pqnn_search
 * 	=========== RICORDARSI DI DEALLOCARE LE COSE OVUNQUE
 */
void pqnn_search(params* input) {
	
	if(input->exaustive==0){
	
		//per ogni punto per query set
		int i,i_w,ind,result,sjump,sbucket;
		c_x=alloc_vector(input->m);
		for(i=0;i< input->nq;i++){
			
			//printf("calcolo dei w centroidi più vicini alla query X = %d\n",i);
			//printX(x_query,i,input->d);
			//calcolo dei w centroidi più vicini a x
			//cerco di passarlgi solo il punto x in modo che i metodi possono preoccuparsi solo di
			//ciclare su 128 dimensioni in quanto punto singolo. (sulle dimensioni in generale)
			//int* label_w =  w_near_centroids(&x_query[i*input->d],Cc,input->kc,input->d,input->w);
			int* label_w = w_near_centroids(&x_query[i*input->d],Cc,input->kc,input->d,input->w);

			//test per il nuovo metodo w_near_centroids
			//printVector(label_w,input->w);
			//printVector(label_w_2,input->w);


			// per ogni centroide vicino appiclo la ricerca
			
			//calcolo tutti i residui r(x) con i centroidi in w
			//printf("calcolo dei residui r(x) con tutti i centroidi w\n");
			float* res_x= residuals_x(&x_query[i*input->d],Cc,label_w,input->w,input->d);

			//funzione NASM per il calcolo dei residui


			//da testare meglio per vedere se
			//effettivamente funziona

			//allocazioni per ottenere un MaxHeap che interagisca con il metodo max_heap
			int * k_nn=alloc_vector(input->knn);
			float * result_dist=alloc_matrix(input->knn,1);
			float tmp,nn_dis = FLT_MAX;//DBL_MAX;
			int C_i;
			int* L_i;
			//printDsQs(res_x,NULL,input->w,input->d,0);
			for(i_w = 0 ; i_w < input->w ; i_w++){
				

				if(input->symmetric==1)
				{
					//printf("SDC: scorrimento della Inverted List: %d\n",i_w);
					for(int j=0;j<input->m;j++){
						//uj_x = Uj( &x_query[i*input->d], j, input->m,1,input->d);
						uj_x = Uj( &res_x[i_w*input->d], j, input->m,1,input->d);
						c_x[j] = centX(&Cp[j*input->sub*input->k], uj_x, input->k, input->d/input->m);
					}	
					dealloc_matrix(uj_x); //capire se è necessario perchè sembra che perda molto tempo
					//centroide più vicino associato al punto
					C_i= label_w[i_w];
					L_i = IL[C_i];

					//variabile usata per non accedere continuamente in jump[C_i]*nodo
					//sjump=jump[C_i]*nodo;
					sbucket=bucket[C_i];
					
					
					

					//calcolo tutte le distanze tra res(x) e le Cji presenti nella inverted List
					for( ind = 0; ind < sbucket; ind++)
					{	
						
						//printVector(&L_i[ind][1],input->m);
						tmp = NE_sdc(c_x,stored_distance, input->m, &L_i[ind*nodo +1],input->k);
						if(tmp < nn_dis){
							
							nn_dis=max_heap(k_nn,result_dist,L_i[ind*nodo],tmp,nn_dis,input->knn,false);
							//printf("\n nn_dis = %f  per il punto y = %d\n\n",nn_dis,L_i[ind][0]);
							//nn_dis = tmp;
							//result = L_i[ind][0];
						}
					}
				}
				else if(input->symmetric==0)
				{
						
					//calcolare tutte le distanze d(Uj(r(x)),Cji)^2
					//per ogni sotto quantizzatore j e per ogni centroide Cji

					//printf("ADC: scorrimento della Inverted List: %d\n",i_w);
					stored_distance=pre_adc(&res_x[i_w*input->d],Cp,input->d,input->m,input->k);

					//adesso devo entrare nell'inverted list con il centroide w' in questione e calcolare 
					//la distanza con tutte le y che fanno parte della lista (che poi servirà avere in centroide
					// collegato a quella determinata y)

					//con questo ottendo la i-esima inverted List 
					C_i= label_w[i_w];
					L_i = IL[C_i];
					//sjump=jump[C_i]*nodo;
					sbucket=bucket[C_i];


					//calcolo tutte le distanze tra res(x) e le Cji presenti nella inverted List
					for( ind = 0; ind < sbucket; ind++)
					{
						tmp = NE_adc(stored_distance, input->k, input->m, &L_i[ind*nodo +1]);
						if(tmp < nn_dis){
							nn_dis=max_heap(k_nn,result_dist, L_i[ind*nodo],tmp,nn_dis,input->knn,false);
							//printf("\n nn_dis = %f  per il punto y = %d\n\n",nn_dis,L_i[ind][0]);
							//nn_dis = tmp;
							//result = L_i[ind][0];
						}
					}
					
				}
				
				
			}
			//salvarsi i k valori che minimizzano la distanza(con il MAXHEAP) non ordinato
			//printf("per il punto x in posizione %d, il nn è la y in posizione %d \n", i, result);
			for(int j = 0; j < input->knn; j++)
			{
				input->ANN[i*input->knn+j]=k_nn[j];
			}

			//pulizia del max_heap (sarebbe buono capire se conviene deallocarli solo alla fine 
			//oppure ogni volta che dobbiamo "resettarli")
			dealloc_vector(k_nn);
			dealloc_matrix(result_dist);
			c_max_heap=0;
			pre_max_heap=0;
			
			//printf("\n x=%d 	y=%d	dist=%f\n",i,result,nn_dis);
		}
	}


	/* Variante senza max_heap
	
	if(input->exaustive==1 && input->symmetric==1){
		float tmp,nn_dis;
		c_x=alloc_vector(input->m);
		int x,y,index,k;
		for(x=0; x<input->nq;x++){
			bool* r = (bool*) get_block(sizeof(bool),input->n);			
			for(int j=0;j<input->m;j++){
				uj_x = Uj( &input->qs[x*input->d], j, input->m,1,input->d);
				c_x[j] = centX(centroids[j], uj_x, input->k, input->d/input->m);
			}	
			dealloc_matrix(uj_x);
			for(k=0; k<input->knn; k++){
				nn_dis = DBL_MAX;
				for(y=0; y< input->n; y++){
					if(r[y]!=true){
						tmp = sdc(c_x,stored_distance, y, input->m, pq, input->k);
						if(tmp < nn_dis){
							nn_dis = tmp;
							index = y;
						}
					}
				}//for y
				input->ANN[x*input->knn+k]=index;
			//	printf("\n x = %d	y = %d	dist = %f",x,index,nn_dis);
				r[index]=true;
			}//for k
		}
	}*/
	if(input->exaustive==1 && input->symmetric==1){
		float tmp,nn_dis;
		c_x=alloc_vector(input->m);
		int x,y,k;
		for(x=0; x<input->nq;x++){	
			//clock_t t11 = clock();
			int* k_nn = alloc_vector(input->knn);	
			float* result_dist=alloc_matrix(input->knn,1);

			for(int j=0;j<input->m;j++){
				uj_x = Uj( &input->qs[x*input->d], j, input->m,1,input->d);
				c_x[j] = centX(&centroids[j*input->k*input->d /input->m], uj_x, input->k, input->d/input->m);
			}	
			dealloc_matrix(uj_x);

			nn_dis = FLT_MAX;//DBL_MAX;
			for(y=0; y< input->n; y++){
				tmp = sdc(c_x,stored_distance, y, input->m, input->n, pq, input->k);
				if(tmp < nn_dis){
					nn_dis = max_heap(k_nn,result_dist,y,tmp,nn_dis,input->knn,false);
					}
			}//for y
			
			/**
			 * stampa risultati con relative distanze
			*/
		/*
			printf("query %d -->",x);
			for(int i=0; i<input->knn; i++){
				printf(" %d	[dist = %f]	",k_nn[i],result_dist[i]);
			}*/

			/**
			 * stampa distanza reale tra x e y
			*/
			/*printf("\nquery %d -->",x);
			for(int i=0; i<input->knn; i++){
				printf(" %d	[dist = %f]	",k_nn[i],dist(&input->qs[x*input->d],&input->ds[k_nn[i]],input->d));
			}*/

		
			for(int i=0; i<input->knn; i++){
				input->ANN[x*input->knn+i]=k_nn[i];
			}
			c_max_heap=0;
			pre_max_heap=0;
			//t11 = clock() - t11;
			//printf("\n tempo di calcolo per x=%d = %.6f secs\n",x, ((float)t11)/CLOCKS_PER_SEC);
		}
		
	}
	if(input->exaustive==1 && input->symmetric==0){
		float tmp,nn_dis;
		int x,y,k;
		for(x=0; x<input->nq;x++){
			//clock_t t11 = clock();
			int* k_nn = alloc_vector(input->knn);
			float* result_dist=alloc_matrix(input->knn,1);

			stored_distance=pre_adc(&input->qs[x*input->d],centroids,input->d,input->m,input->k);
			nn_dis = FLT_MAX;
			//nn_dis = DBL_MAX;
			for(y=0; y< input->n; y++){
					tmp = adc(stored_distance, y, input->k, input->m, input->n, pq);
					if(tmp < nn_dis){
						nn_dis = max_heap(k_nn,result_dist,y,tmp,nn_dis,input->knn,false);
					}
			}
				/**
			 	* stampa risultati con relative distanze
				*/
				/*printf("query %d -->",x);
				for(int i=0; i<input->knn; i++){
					printf(" %d	[dist = %f]	",k_nn[i],result_dist[i]);
				}
				printf("\nquery %d -->",x);
			for(int i=0; i<input->knn; i++){
				printf(" %d	[dist = %f]	",k_nn[i],dist(&input->qs[x*input->d],&input->ds[k_nn[i]],input->d));
			}
			printf("\n \n");
			*/
				for(int i=0; i<input->knn; i++){
					input->ANN[x*input->knn+i]=k_nn[i];
				}
				c_max_heap=0;
				pre_max_heap = 0;	
				//t11 = clock() - t11;
			//printf("\n tempo di calcolo per x=%d = %.6f secs\n",x, ((float)t11)/CLOCKS_PER_SEC);
		}
		
	}

	/*printf("\nACCESSI SDC : %d\n",accesso);
	printf("\nACCESSI ADC : %d\n",accesso_1);
	printf("\nACCESSI nel max_heap : %d\n",accesso_2);
	*/


	/* variante senza max_heap
	if(input->exaustive==1 && input->symmetric==0){
		float tmp,nn_dis;
		int x,y,index,k;
		for(x=0; x<input->nq;x++){
			bool* r = (bool*) get_block(sizeof(bool),input->n);			
			//printf("Calcolo le distanze (ASIMMETRICO) tra X e Cji\n");
			stored_distance=pre_adc(&input->qs[x*input->d],centroids,input->d,input->m,input->k);
			for(k=0; k<input->knn; k++){
				nn_dis = DBL_MAX;
				for(y=0; y< input->n; y++){
					if(r[y]!=true){
						tmp = adc(stored_distance, y, input->m, pq);
						if(x==999){
							printf("y = %d 	dist = %f \n",y,tmp);
						}
						if(tmp < nn_dis){
							nn_dis = tmp;
							index = y;
						}
					}
				}
				input->ANN[x*input->knn+k]=index;
			//	printf("\n x = %d	y = %d	dist = %f",x,index,nn_dis);
				r[index]=true;	
			}
		}
	}*/
	
    // -------------------------------------------------
    // Codificare qui l'algoritmo di interrogazione
    // -------------------------------------------------
    
    //pqnn32_search(input); // Chiamata funzione assembly

	// Restituisce il risultato come una matrice di nq * knn
	// identificatori associati agli ANN approssimati delle nq query.
	// La matrice è memorizzata per righe
    // -------------------------------------------------

}


int main(int argc, char** argv) {

	
	
	char fname[256];
	int i, j;
	
	//
	// Imposta i valori di default dei parametri
	//

	params* input = malloc(sizeof(params));

	input->filename = NULL;
	input->exaustive = 1;
	input->symmetric = 1;
	input->knn = 4;
	input->m = 8;
	input->k = 32;
	//input->kc = 8192;
	input->kc = 32;
	//input->w = 16;
	input->w=8;
	input->eps = 0.01;
	input->tmin = 10;
	input->tmax = 100;
	input->silent = 0;
	input->display = 0;

	//
	// Legge i valori dei parametri da riga comandi
	//

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
	
	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

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
	
	//
	// Legge il data set ed il query set
	//
	
	if (input->filename == NULL || strlen(input->filename) == 0) {
		printf("Missing input file name!\n");
		exit(1);
	}
	
	sprintf(fname, "%s.ds", input->filename);
	input->ds = load_data(fname, &input->n, &input->d);
	input->sub=input->d/input->m;
	//input->n = input->n/2 + 2;

	input->nr = input->n/20;

	sprintf(fname, "%s.qs", input->filename);
	input->qs = load_data(fname, &input->nq, &input->d);
	//input->nq=input->nq/2 + 2;

	//creazione di una matrice temporanea che ospita un sottogruppo di dimensioni del dataset (n*sub dimensionale)
	
	/*MATRIX tmp = Uj(input->ds,0,input->m,input->n,input->d);
	printEq(input->ds,tmp,input->n,input->d,input->n,input->sub);
	printf("\n QUERY SET");
	*/

	//
	// Visualizza il valore dei parametri
	//
	
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
	
	//
	// Costruisce i quantizzatori
	//
	
	clock_t t = clock();
	pqnn_index(input);
	t = clock() - t;
	
	if (!input->silent)
		printf("\nIndexing time = %.3f secs\n", ((float)t)/CLOCKS_PER_SEC);
	else
		printf("%.3f\n", ((float)t)/CLOCKS_PER_SEC);

	//
	// Determina gli ANN
	//
	
	input->ANN = calloc(input->nq*input->knn,sizeof(int));

	clock_t t_1 = clock();
	pqnn_search(input);
	t_1 = clock() - t_1;
	
	if (!input->silent)
		printf("\nSearching time = %.3f secs\n", ((float)t_1)/CLOCKS_PER_SEC);
	else
		printf("%.3f\n", ((float)t_1)/CLOCKS_PER_SEC);
	
	//
	// Salva gli ANN
	//
	
 	if (input->ANN != NULL)
 	{
 		//if (!input->silent && input->display) { se si discommenta ritorna "opzionale"
 			printf("\nANN:\n");
 			for (i = 0; i < input->nq; i++) {
				printf("query #%d:", i);
				for (j = 0; j < input->knn; j++)
					printf("  %d  ", input->ANN[i*input->knn+j]);
				printf("\n");
 			}
 		//}
 		save_ANN(input->filename, input->ANN, input->nq, input->knn);
 	}
	
	if (!input->silent){
		printf("\nDone.\n");
	}

	printf("\nIndexing time = %.3f secs\n", ((float)t)/CLOCKS_PER_SEC);
	printf("\nSearching time = %.3f secs\n", ((float)t_1)/CLOCKS_PER_SEC);
	//printDsQs(input->ds,NULL,input->n,input->d,0);
	return 0;
}
