#include "matpr.h"

/* double precision matrix print function: */
void matrix_print(char *name, double *a, int *m, int *n, int *lda, int *level){

    int i,j;
    FILE *file;

    switch (*level) {
    case 2: 
	/* stdout output: */
	for (i=0; i<*m; i++){
	    for (j=0; j<*n; j++){
		printf("%f8.4 ", a[i+j*(*lda)]);
	    }
	    printf("\n");
	}
    default:
    case 1:
	/* create octave matrix file: */
	file=fopen(name,"w");
	if(!file){
	    printf("error opening file");
	    break;
	}
	fprintf(file, "# Created by matpr\n");
	fprintf(file, "# name: %s\n", name);
	fprintf(file, "# type: matrix\n");
	fprintf(file, "# rows: %i\n", *m);
	fprintf(file, "# columns: %i\n", *n);
	fflush(file);
	for (i=0; i<*m; i++){
	    for (j=0; j<*n; j++){
		fprintf(file," %f ", a[i+j*(*lda)]);
	    }
	    fprintf(file,"\n");
	    }	
	fclose(file);
	
    }
}

/* Fortran wrapper: */
void F77_CALL(matpr) (char *name, double *a, int *m, int *n, int *lda, int *level){
    matrix_print(name, a, m, n, lda, level);
}
/* integer matrix print function: */
void imatrix_print(char *name, int *a, int *m, int *n, int *lda, int *level){

    int i,j;
    FILE *file;

    switch (*level) {
    case 2: 
	/* stdout output: */
	for (i=0; i<*m; i++){
	    for (j=0; j<*n; j++){
		printf("%i ", a[i+j*(*lda)]);
	    }
	    printf("\n");
	}
    default:
    case 1:
	/* create octave matrix file: */
	file=fopen(name,"w");
	if(!file){
	    printf("error opening file");
	    break;
	}
	fprintf(file, "# Created by imatpr\n");
	fprintf(file, "# name: %s\n", name);
	fprintf(file, "# type: matrix\n");
	fprintf(file, "# rows: %i\n", *m);
	fprintf(file, "# columns: %i\n", *n);
	fflush(file);
	for (i=0; i<*m; i++){
	    for (j=0; j<*n; j++){
		fprintf(file," %i ", a[i+j*(*lda)]);
	    }
	    fprintf(file,"\n");
	    }	
	fclose(file);
	
    }
}

/* Fortran wrapper: */
void F77_CALL(imatpr) (char *name, int *a, int *m, int *n, int *lda, int *level){
    imatrix_print(name, a, m, n, lda, level);
}
