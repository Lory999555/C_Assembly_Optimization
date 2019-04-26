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
//#include "kmeans2.c"

#define 	MATRIX		double*
#define	    VECTOR		double*


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
	//
	int* ANN; 
	//
	// Inserire qui i campi necessari a memorizzare i Quantizzatori
	//
	// ...
	// ...
	// ...
	//
} params;


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
	return (MATRIX) get_block(sizeof(double),rows*cols);
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
	status = fread(data, sizeof(double), rows*cols, fp);
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
			fprintf(fp, " %f  ", ANN[i*knn+j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
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

//metodo che prende un sottogruppo (sub dimensionale) del data set
// j serve per prendere il j-esimo gruppetto di dimensioni j=2 equivale ad U2
MATRIX Uj(MATRIX ds, int j,int m,int n,int d){
	int i,sub=d/m;
	int k,c;
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

double dist(double * x,double * y, int k){
	double distance = 0;
	for (int i=0; i<k;i++){
	    distance += pow(x[i] - y[i], 2);
	}
	return distance;
}

MATRIX randCentroid(MATRIX ds,int n,int d,int k){
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
}


//kmeans modificato in modo da prendere due "MATRIX" e usando l'alloc del prof con l'allineamento.
int * k_means(double * data, int n, int d, int k, double t, double * centroids,int t_min,int t_max) {
	
	/* output cluster label for each data point */
	int * labels = alloc_vector(n);

	double min_distance;
	double distance;
	double offset;
	int iter=0;
	int h, i, j; /* loop counters, of course */
	
	/* size of each cluster */
	int* counts = alloc_vector(k);
	double old_error, error = DBL_MAX; /* sum of squared euclidean distance */
	
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
			printf("c[%d][%d] = data[%d][%d] ------>  %f  ||||  %f \n",i,j,h,j,c[i*d+j],data[h*d+j]);
		}
	}
	printf("main loop\n");
	/****
	 ** main loop */

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
			min_distance = DBL_MAX;
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
					printf("SI ---> c[%d][%d] =  %f \n",i,j,c[i*d+j]);
				}else
				{
					c[i*d+j] =c1[i*d+j];
					printf("NO ---> c[%d][%d] =  %f \n",i,j,c[i*d+j]);
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

	return labels;
}




extern void pqnn32_index(params* input);
extern int* pqnn32_search(params* input);



/*
 *	pqnn_index
 * 	==========
 */
void pqnn_index(params* input) {

	//printDsQs(input->ds,NULL,input->n,input->d,0);
	if(input->exaustive == 0){
		/*
		era per provare il kmeans1.c
		MATRIX Cc = randCentroid(input->ds,input->n,input->d,input->kc);
		printDsQs(Cc,NULL,input->kc,input->d,0);
		int* Cc_index = (int*)get_block(sizeof(int),input->n);
		kmeans(input->d,input->ds,input->n,input->kc,Cc,Cc_index);
		printCentroids(Cc,Cc_index,input->n,input->d,input->kc);
		*/
		MATRIX Cc= alloc_matrix(input->kc,input->d);
		int * Cc_index=k_means(input->ds,input->n,input->d,input->kc,input->eps,Cc,input->tmin,input->tmax);
		printCentroids(Cc,Cc_index,input->n,input->d,input->kc);
	}
    // -------------------------------------------------
    // Codificare qui l'algoritmo di indicizzazione
    // -------------------------------------------------
    
    pqnn32_index(input); // Chiamata funzione assembly

    // -------------------------------------------------

}


/*
 *	pqnn_search
 * 	===========
 */
void pqnn_search(params* input) {
	
    // -------------------------------------------------
    // Codificare qui l'algoritmo di interrogazione
    // -------------------------------------------------
    
    pqnn32_search(input); // Chiamata funzione assembly

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
	input->knn = 1;
	input->m = 8;
	input->k = 256;
	//input->kc = 8192;
	input->kc = 16;
	input->w = 16;
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

	input->nr = input->n/20;

	sprintf(fname, "%s.qs", input->filename);
	input->qs = load_data(fname, &input->nq, &input->d);

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

	t = clock();
	pqnn_search(input);
	t = clock() - t;
	
	if (!input->silent)
		printf("\nSearching time = %.3f secs\n", ((float)t)/CLOCKS_PER_SEC);
	else
		printf("%.3f\n", ((float)t)/CLOCKS_PER_SEC);
	
	//
	// Salva gli ANN
	//
	
 	if (input->ANN != NULL)
 	{
 		if (!input->silent && input->display) {
 			printf("\nANN:\n");
 			/*for (i = 0; i < input->nq; i++) {
				printf("query #%d:", i);
				for (j = 0; j < input->knn; j++)
					printf("  %f ", input->ANN[i*input->knn+j]);
				printf("\n");
 			}*/
 		}
 		save_ANN(input->filename, input->ANN, input->nq, input->knn);
	}
	
	if (!input->silent)
		printf("\nDone.\n");

	return 0;
}
