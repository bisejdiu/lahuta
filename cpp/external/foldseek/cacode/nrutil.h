// Public domain code from Yi-Kuo Yu & Stephen Altschul, NCBI

float *vector();
float **matrix();
float **convert_matrix();
double *dvector();
double *dvector(int nl, int nh);
double **dmatrix();
double **dmatrix(int nrl, int nrh, int ncl, int nch);
int *ivector();
int *ivector(int nl, int nh);
int **imatrix();
float **submatrix();
void free_vector();
void free_dvector();
void free_dvector(double *v, int nl, int nh);
void free_ivector();
void free_ivector(int *v, int nl, int nh);
void free_matrix();
void free_dmatrix();
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix();
void free_submatrix();
void free_convert_matrix();
void nrerror();
void nrerror(const char *error_text);
