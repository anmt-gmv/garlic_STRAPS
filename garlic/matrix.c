#include "matrix.h"

/* Returns the vector index referred to a set (row,col) in the matrix */
static int get_sym_index(int row, int col)
{
	if (col <= row)
	{
		return ((row * (row + 1)) >> 1) + col;
	}
	else
	{
		return ((col * (col + 1)) >> 1) + row;
	}
}

/* Returns the value of the element in matrix in (row,col) */
double get_sym_element(const float *matrix, int row, int col)
{
	int index = get_sym_index(row, col);
	return matrix[index];
}

/* Set the value of the element in matrix in (row,col) */
void set_sym_element(float *matrix, int row, int col, float value)
{
	int index = get_sym_index(row, col);
	matrix[index] = value;
}

/* Adds a value to the element in matrix in (row,col) */
void update_sym_element(float *matrix, int row, int col, float value)
{
	int index      = get_sym_index(row, col);
	matrix[index] += value;
}

/* This function updates the covariance matrix in a Kalman Filter P+ = (I-KH)*P- */
void update_matrix(int dim, const float *K, const float *PH, float *P)
{
	int i, j;
	int index = 0;
	for (i = 0; i < dim; i++)
	{
		float aux = -K[i];
		for (j = 0; j <= i; j++)
		{
			P[index] += aux * PH[j];
			index++;
		}
	}
}

/* This function performs the multiplication between a symmetric matrix A and a vector X: Y = AX */
void matSymVecMul_f(float *Y, const float *A, const float *X, int n)
{
	// Loop variables.
	int i, j;

	// Initialization.
	for(i = 0; i < n; i++) { Y[i] = 0.0; }

	// Auxiliar index.
	int index = 0;
	// Calculate the product.
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < i; j++)
		{
			Y[i] += A[index] * X[j];
			Y[j] += A[index] * X[i];
			index++;
		}
		Y[i] += A[index] * X[i];
		index++;
	}
}

/* This function initializes the sparse matrix */
void reset_sparse_matrix(unsigned int size, sparse_matrix_t *matrix)
{
	int k;
	for (k = 0; k < size; k++)
	{
		matrix->num_elem_per_row[k] = 0;
		matrix->diagonal[k]         = 1.0;
	}
}

/* This function adds an element in the sparse matrix, updating the sparse structure properly */
void add_element_to_sparse_matrix(int row, int col, float value, sparse_matrix_t *matrix)
{
	if (row == col)
	{
		matrix->diagonal[row] += value;
	}
	else
	{
		int index = matrix->num_elem_per_row[row];

		matrix->num_elem_per_row[row]++;

		matrix->values     [row][index] = value;
		matrix->col_number [row][index] = col;
	}
}

/* This function is used to propagate the covariance matrix in the Kalman Filter. The
 * current implementation shall only be invoked when working with the 18-states KF,
 * since it has been optimized to work faster than the nominal and generic function.
 * In this case, the sparse matrix is the transition matrix F, P(k) = F*P(k-1)*F'
 */
#if USE_FAST_MECHANIZATION == 1
void update_matrix_sparse_fast(const sparse_matrix_t *matrix, float *P_mat)
{
	float cov_mat [KNSTATES][KNSTATES];
	float PF_mat  [KNSTATES][KNSTATES-6];

	int i, j;

	int index = 0;
	for (i = 0; i < KNSTATES; i++)
	{
		for (j = 0; j <= i; j++)
		{
			cov_mat[i][j] = P_mat[index];
			cov_mat[j][i] = P_mat[index];
			index++;
		}
	}

	for (i = 0; i < KNSTATES; i++)
	{
#if KAC == 0
		PF_mat[i][0]   = cov_mat[i][3]  * matrix->values[0][0];
		PF_mat[i][1]   = cov_mat[i][4]  * matrix->values[1][0];
		PF_mat[i][2]   = cov_mat[i][5]  * matrix->values[2][0];
		PF_mat[i][3]   = cov_mat[i][0]  * matrix->values[3][0];
		PF_mat[i][3]  += cov_mat[i][4]  * matrix->values[3][1];
		PF_mat[i][3]  += cov_mat[i][10] * matrix->values[3][2];
		PF_mat[i][3]  += cov_mat[i][11] * matrix->values[3][3];
		PF_mat[i][3]  += cov_mat[i][12] * matrix->values[3][4];
		PF_mat[i][3]  += cov_mat[i][13] * matrix->values[3][5];
		PF_mat[i][3]  += cov_mat[i][14] * matrix->values[3][6];
		PF_mat[i][4]   = cov_mat[i][1]  * matrix->values[4][0];
		PF_mat[i][4]  += cov_mat[i][3]  * matrix->values[4][1];
		PF_mat[i][4]  += cov_mat[i][9]  * matrix->values[4][2];
		PF_mat[i][4]  += cov_mat[i][11] * matrix->values[4][3];
		PF_mat[i][4]  += cov_mat[i][12] * matrix->values[4][4];
		PF_mat[i][4]  += cov_mat[i][13] * matrix->values[4][5];
		PF_mat[i][4]  += cov_mat[i][14] * matrix->values[4][6];
		PF_mat[i][5]   = cov_mat[i][9]  * matrix->values[5][0];
		PF_mat[i][5]  += cov_mat[i][10] * matrix->values[5][1];
		PF_mat[i][5]  += cov_mat[i][12] * matrix->values[5][2];
		PF_mat[i][5]  += cov_mat[i][13] * matrix->values[5][3];
		PF_mat[i][5]  += cov_mat[i][14] * matrix->values[5][4];
		PF_mat[i][6]   = cov_mat[i][8]  * matrix->values[6][0];
		PF_mat[i][7]   = 0.0;
		PF_mat[i][8]   = 0.0;
		PF_mat[i][9]   = cov_mat[i][10] * matrix->values[9][0];
		PF_mat[i][9]  += cov_mat[i][15] * matrix->values[9][1];
		PF_mat[i][9]  += cov_mat[i][16] * matrix->values[9][2];
		PF_mat[i][9]  += cov_mat[i][17] * matrix->values[9][3];
		PF_mat[i][10]  = cov_mat[i][9]  * matrix->values[10][0];
		PF_mat[i][10] += cov_mat[i][15] * matrix->values[10][1];
		PF_mat[i][10] += cov_mat[i][16] * matrix->values[10][2];
		PF_mat[i][10] += cov_mat[i][17] * matrix->values[10][3];
		PF_mat[i][11]  = cov_mat[i][15] * matrix->values[11][0];
		PF_mat[i][11] += cov_mat[i][16] * matrix->values[11][1];
		PF_mat[i][11] += cov_mat[i][17] * matrix->values[11][2];
#else
		PF_mat[i][0]   = cov_mat[3][i]  * matrix->values[0][0];
		PF_mat[i][1]   = cov_mat[4][i]  * matrix->values[1][0];
		PF_mat[i][2]   = cov_mat[5][i]  * matrix->values[2][0];
		PF_mat[i][3]   = cov_mat[0][i]  * matrix->values[3][0];
		PF_mat[i][3]  += cov_mat[4][i]  * matrix->values[3][1];
		PF_mat[i][3]  += cov_mat[11][i] * matrix->values[3][2];
		PF_mat[i][3]  += cov_mat[12][i] * matrix->values[3][3];
		PF_mat[i][3]  += cov_mat[13][i] * matrix->values[3][4];
		PF_mat[i][3]  += cov_mat[14][i] * matrix->values[3][5];
		PF_mat[i][3]  += cov_mat[15][i] * matrix->values[3][6];
		PF_mat[i][4]   = cov_mat[1][i]  * matrix->values[4][0];
		PF_mat[i][4]  += cov_mat[3][i]  * matrix->values[4][1];
		PF_mat[i][4]  += cov_mat[10][i] * matrix->values[4][2];
		PF_mat[i][4]  += cov_mat[12][i] * matrix->values[4][3];
		PF_mat[i][4]  += cov_mat[13][i] * matrix->values[4][4];
		PF_mat[i][4]  += cov_mat[14][i] * matrix->values[4][5];
		PF_mat[i][4]  += cov_mat[15][i] * matrix->values[4][6];
		PF_mat[i][5]   = cov_mat[10][i] * matrix->values[5][0];
		PF_mat[i][5]  += cov_mat[11][i] * matrix->values[5][1];
		PF_mat[i][5]  += cov_mat[13][i] * matrix->values[5][2];
		PF_mat[i][5]  += cov_mat[14][i] * matrix->values[5][3];
		PF_mat[i][5]  += cov_mat[15][i] * matrix->values[5][4];
		PF_mat[i][6]   = cov_mat[9][i]  * matrix->values[6][0];
		PF_mat[i][7]   = 0;
		PF_mat[i][8]   = 0;
		PF_mat[i][9]   = 0;
		PF_mat[i][10]  = cov_mat[11][i] * matrix->values[10][0];
		PF_mat[i][10] += cov_mat[16][i] * matrix->values[10][1];
		PF_mat[i][10] += cov_mat[17][i] * matrix->values[10][2];
		PF_mat[i][10] += cov_mat[18][i] * matrix->values[10][3];
		PF_mat[i][11]  = cov_mat[10][i] * matrix->values[11][0];
		PF_mat[i][11] += cov_mat[16][i] * matrix->values[11][1];
		PF_mat[i][11] += cov_mat[17][i] * matrix->values[11][2];
		PF_mat[i][11] += cov_mat[18][i] * matrix->values[11][3];
		PF_mat[i][12]  = cov_mat[16][i] * matrix->values[12][0];
		PF_mat[i][12] += cov_mat[17][i] * matrix->values[12][1];
		PF_mat[i][12] += cov_mat[18][i] * matrix->values[12][2];
#endif
	}

	float *p_ptr = &P_mat[0];
	for (i = 0; i < KNSTATES-6; i++)
	{
		for (j = 0; j <= i; j++)
		{
			p_ptr[j] += PF_mat[j][i];
			p_ptr[j] += PF_mat[i][j];
		}

		p_ptr += i+1;
	}

	for (; i < KNSTATES; i++)
	{
		for (j = 0; j < KNSTATES-6; j++)
		{
			p_ptr[j] *= matrix->diagonal[i];
			p_ptr[j] += PF_mat[i][j];
		}

		for (; j <= i; j++)
		{
			p_ptr[j] *= matrix->diagonal[i] * matrix->diagonal[j];
		}

		p_ptr += i+1;
	}
#if KAC == 0
	P_mat[0]  += matrix->values[0][0] * PF_mat[3][0];
	P_mat[1]  += matrix->values[1][0] * PF_mat[4][0];
	P_mat[2]  += matrix->values[1][0] * PF_mat[4][1];
	P_mat[3]  += matrix->values[2][0] * PF_mat[5][0];
	P_mat[4]  += matrix->values[2][0] * PF_mat[5][1];
	P_mat[5]  += matrix->values[2][0] * PF_mat[5][2];
	P_mat[6]  += matrix->values[3][0] * PF_mat[0][0];
	P_mat[7]  += matrix->values[3][0] * PF_mat[0][1];
	P_mat[8]  += matrix->values[3][0] * PF_mat[0][2];
	P_mat[9]  += matrix->values[3][0] * PF_mat[0][3];
	P_mat[6]  += matrix->values[3][1] * PF_mat[4][0];
	P_mat[7]  += matrix->values[3][1] * PF_mat[4][1];
	P_mat[8]  += matrix->values[3][1] * PF_mat[4][2];
	P_mat[9]  += matrix->values[3][1] * PF_mat[4][3];
	P_mat[6]  += matrix->values[3][2] * PF_mat[10][0];
	P_mat[7]  += matrix->values[3][2] * PF_mat[10][1];
	P_mat[8]  += matrix->values[3][2] * PF_mat[10][2];
	P_mat[9]  += matrix->values[3][2] * PF_mat[10][3];
	P_mat[6]  += matrix->values[3][3] * PF_mat[11][0];
	P_mat[7]  += matrix->values[3][3] * PF_mat[11][1];
	P_mat[8]  += matrix->values[3][3] * PF_mat[11][2];
	P_mat[9]  += matrix->values[3][3] * PF_mat[11][3];
	P_mat[6]  += matrix->values[3][4] * PF_mat[12][0];
	P_mat[7]  += matrix->values[3][4] * PF_mat[12][1];
	P_mat[8]  += matrix->values[3][4] * PF_mat[12][2];
	P_mat[9]  += matrix->values[3][4] * PF_mat[12][3];
	P_mat[6]  += matrix->values[3][5] * PF_mat[13][0];
	P_mat[7]  += matrix->values[3][5] * PF_mat[13][1];
	P_mat[8]  += matrix->values[3][5] * PF_mat[13][2];
	P_mat[9]  += matrix->values[3][5] * PF_mat[13][3];
	P_mat[6]  += matrix->values[3][6] * PF_mat[14][0];
	P_mat[7]  += matrix->values[3][6] * PF_mat[14][1];
	P_mat[8]  += matrix->values[3][6] * PF_mat[14][2];
	P_mat[9]  += matrix->values[3][6] * PF_mat[14][3];
	P_mat[10] += matrix->values[4][0] * PF_mat[1][0];
	P_mat[11] += matrix->values[4][0] * PF_mat[1][1];
	P_mat[12] += matrix->values[4][0] * PF_mat[1][2];
	P_mat[13] += matrix->values[4][0] * PF_mat[1][3];
	P_mat[14] += matrix->values[4][0] * PF_mat[1][4];
	P_mat[10] += matrix->values[4][1] * PF_mat[3][0];
	P_mat[11] += matrix->values[4][1] * PF_mat[3][1];
	P_mat[12] += matrix->values[4][1] * PF_mat[3][2];
	P_mat[13] += matrix->values[4][1] * PF_mat[3][3];
	P_mat[14] += matrix->values[4][1] * PF_mat[3][4];
	P_mat[10] += matrix->values[4][2] * PF_mat[9][0];
	P_mat[11] += matrix->values[4][2] * PF_mat[9][1];
	P_mat[12] += matrix->values[4][2] * PF_mat[9][2];
	P_mat[13] += matrix->values[4][2] * PF_mat[9][3];
	P_mat[14] += matrix->values[4][2] * PF_mat[9][4];
	P_mat[10] += matrix->values[4][3] * PF_mat[11][0];
	P_mat[11] += matrix->values[4][3] * PF_mat[11][1];
	P_mat[12] += matrix->values[4][3] * PF_mat[11][2];
	P_mat[13] += matrix->values[4][3] * PF_mat[11][3];
	P_mat[14] += matrix->values[4][3] * PF_mat[11][4];
	P_mat[10] += matrix->values[4][4] * PF_mat[12][0];
	P_mat[11] += matrix->values[4][4] * PF_mat[12][1];
	P_mat[12] += matrix->values[4][4] * PF_mat[12][2];
	P_mat[13] += matrix->values[4][4] * PF_mat[12][3];
	P_mat[14] += matrix->values[4][4] * PF_mat[12][4];
	P_mat[10] += matrix->values[4][5] * PF_mat[13][0];
	P_mat[11] += matrix->values[4][5] * PF_mat[13][1];
	P_mat[12] += matrix->values[4][5] * PF_mat[13][2];
	P_mat[13] += matrix->values[4][5] * PF_mat[13][3];
	P_mat[14] += matrix->values[4][5] * PF_mat[13][4];
	P_mat[10] += matrix->values[4][6] * PF_mat[14][0];
	P_mat[11] += matrix->values[4][6] * PF_mat[14][1];
	P_mat[12] += matrix->values[4][6] * PF_mat[14][2];
	P_mat[13] += matrix->values[4][6] * PF_mat[14][3];
	P_mat[14] += matrix->values[4][6] * PF_mat[14][4];
	P_mat[15] += matrix->values[5][0] * PF_mat[9][0];
	P_mat[16] += matrix->values[5][0] * PF_mat[9][1];
	P_mat[17] += matrix->values[5][0] * PF_mat[9][2];
	P_mat[18] += matrix->values[5][0] * PF_mat[9][3];
	P_mat[19] += matrix->values[5][0] * PF_mat[9][4];
	P_mat[20] += matrix->values[5][0] * PF_mat[9][5];
	P_mat[15] += matrix->values[5][1] * PF_mat[10][0];
	P_mat[16] += matrix->values[5][1] * PF_mat[10][1];
	P_mat[17] += matrix->values[5][1] * PF_mat[10][2];
	P_mat[18] += matrix->values[5][1] * PF_mat[10][3];
	P_mat[19] += matrix->values[5][1] * PF_mat[10][4];
	P_mat[20] += matrix->values[5][1] * PF_mat[10][5];
	P_mat[15] += matrix->values[5][2] * PF_mat[12][0];
	P_mat[16] += matrix->values[5][2] * PF_mat[12][1];
	P_mat[17] += matrix->values[5][2] * PF_mat[12][2];
	P_mat[18] += matrix->values[5][2] * PF_mat[12][3];
	P_mat[19] += matrix->values[5][2] * PF_mat[12][4];
	P_mat[20] += matrix->values[5][2] * PF_mat[12][5];
	P_mat[15] += matrix->values[5][3] * PF_mat[13][0];
	P_mat[16] += matrix->values[5][3] * PF_mat[13][1];
	P_mat[17] += matrix->values[5][3] * PF_mat[13][2];
	P_mat[18] += matrix->values[5][3] * PF_mat[13][3];
	P_mat[19] += matrix->values[5][3] * PF_mat[13][4];
	P_mat[20] += matrix->values[5][3] * PF_mat[13][5];
	P_mat[15] += matrix->values[5][4] * PF_mat[14][0];
	P_mat[16] += matrix->values[5][4] * PF_mat[14][1];
	P_mat[17] += matrix->values[5][4] * PF_mat[14][2];
	P_mat[18] += matrix->values[5][4] * PF_mat[14][3];
	P_mat[19] += matrix->values[5][4] * PF_mat[14][4];
	P_mat[20] += matrix->values[5][4] * PF_mat[14][5];
	P_mat[21] += matrix->values[6][0] * PF_mat[8][0];
	P_mat[22] += matrix->values[6][0] * PF_mat[8][1];
	P_mat[23] += matrix->values[6][0] * PF_mat[8][2];
	P_mat[24] += matrix->values[6][0] * PF_mat[8][3];
	P_mat[25] += matrix->values[6][0] * PF_mat[8][4];
	P_mat[26] += matrix->values[6][0] * PF_mat[8][5];
	P_mat[27] += matrix->values[6][0] * PF_mat[8][6];
	P_mat[45] += matrix->values[9][0] * PF_mat[10][0];
	P_mat[46] += matrix->values[9][0] * PF_mat[10][1];
	P_mat[47] += matrix->values[9][0] * PF_mat[10][2];
	P_mat[48] += matrix->values[9][0] * PF_mat[10][3];
	P_mat[49] += matrix->values[9][0] * PF_mat[10][4];
	P_mat[50] += matrix->values[9][0] * PF_mat[10][5];
	P_mat[51] += matrix->values[9][0] * PF_mat[10][6];
	P_mat[54] += matrix->values[9][0] * PF_mat[10][9];
	P_mat[45] += matrix->values[9][1] * PF_mat[15][0];
	P_mat[46] += matrix->values[9][1] * PF_mat[15][1];
	P_mat[47] += matrix->values[9][1] * PF_mat[15][2];
	P_mat[48] += matrix->values[9][1] * PF_mat[15][3];
	P_mat[49] += matrix->values[9][1] * PF_mat[15][4];
	P_mat[50] += matrix->values[9][1] * PF_mat[15][5];
	P_mat[51] += matrix->values[9][1] * PF_mat[15][6];
	P_mat[54] += matrix->values[9][1] * PF_mat[15][9];
	P_mat[45] += matrix->values[9][2] * PF_mat[16][0];
	P_mat[46] += matrix->values[9][2] * PF_mat[16][1];
	P_mat[47] += matrix->values[9][2] * PF_mat[16][2];
	P_mat[48] += matrix->values[9][2] * PF_mat[16][3];
	P_mat[49] += matrix->values[9][2] * PF_mat[16][4];
	P_mat[50] += matrix->values[9][2] * PF_mat[16][5];
	P_mat[51] += matrix->values[9][2] * PF_mat[16][6];
	P_mat[54] += matrix->values[9][2] * PF_mat[16][9];
	P_mat[45] += matrix->values[9][3] * PF_mat[17][0];
	P_mat[46] += matrix->values[9][3] * PF_mat[17][1];
	P_mat[47] += matrix->values[9][3] * PF_mat[17][2];
	P_mat[48] += matrix->values[9][3] * PF_mat[17][3];
	P_mat[49] += matrix->values[9][3] * PF_mat[17][4];
	P_mat[50] += matrix->values[9][3] * PF_mat[17][5];
	P_mat[51] += matrix->values[9][3] * PF_mat[17][6];
	P_mat[54] += matrix->values[9][3] * PF_mat[17][9];
	P_mat[55] += matrix->values[10][0] * PF_mat[9][0];
	P_mat[56] += matrix->values[10][0] * PF_mat[9][1];
	P_mat[57] += matrix->values[10][0] * PF_mat[9][2];
	P_mat[58] += matrix->values[10][0] * PF_mat[9][3];
	P_mat[59] += matrix->values[10][0] * PF_mat[9][4];
	P_mat[60] += matrix->values[10][0] * PF_mat[9][5];
	P_mat[61] += matrix->values[10][0] * PF_mat[9][6];
	P_mat[64] += matrix->values[10][0] * PF_mat[9][9];
	P_mat[65] += matrix->values[10][0] * PF_mat[9][10];
	P_mat[55] += matrix->values[10][1] * PF_mat[15][0];
	P_mat[56] += matrix->values[10][1] * PF_mat[15][1];
	P_mat[57] += matrix->values[10][1] * PF_mat[15][2];
	P_mat[58] += matrix->values[10][1] * PF_mat[15][3];
	P_mat[59] += matrix->values[10][1] * PF_mat[15][4];
	P_mat[60] += matrix->values[10][1] * PF_mat[15][5];
	P_mat[61] += matrix->values[10][1] * PF_mat[15][6];
	P_mat[64] += matrix->values[10][1] * PF_mat[15][9];
	P_mat[65] += matrix->values[10][1] * PF_mat[15][10];
	P_mat[55] += matrix->values[10][2] * PF_mat[16][0];
	P_mat[56] += matrix->values[10][2] * PF_mat[16][1];
	P_mat[57] += matrix->values[10][2] * PF_mat[16][2];
	P_mat[58] += matrix->values[10][2] * PF_mat[16][3];
	P_mat[59] += matrix->values[10][2] * PF_mat[16][4];
	P_mat[60] += matrix->values[10][2] * PF_mat[16][5];
	P_mat[61] += matrix->values[10][2] * PF_mat[16][6];
	P_mat[64] += matrix->values[10][2] * PF_mat[16][9];
	P_mat[65] += matrix->values[10][2] * PF_mat[16][10];
	P_mat[55] += matrix->values[10][3] * PF_mat[17][0];
	P_mat[56] += matrix->values[10][3] * PF_mat[17][1];
	P_mat[57] += matrix->values[10][3] * PF_mat[17][2];
	P_mat[58] += matrix->values[10][3] * PF_mat[17][3];
	P_mat[59] += matrix->values[10][3] * PF_mat[17][4];
	P_mat[60] += matrix->values[10][3] * PF_mat[17][5];
	P_mat[61] += matrix->values[10][3] * PF_mat[17][6];
	P_mat[64] += matrix->values[10][3] * PF_mat[17][9];
	P_mat[65] += matrix->values[10][3] * PF_mat[17][10];
	P_mat[66] += matrix->values[11][0] * PF_mat[15][0];
	P_mat[67] += matrix->values[11][0] * PF_mat[15][1];
	P_mat[68] += matrix->values[11][0] * PF_mat[15][2];
	P_mat[69] += matrix->values[11][0] * PF_mat[15][3];
	P_mat[70] += matrix->values[11][0] * PF_mat[15][4];
	P_mat[71] += matrix->values[11][0] * PF_mat[15][5];
	P_mat[72] += matrix->values[11][0] * PF_mat[15][6];
	P_mat[75] += matrix->values[11][0] * PF_mat[15][9];
	P_mat[76] += matrix->values[11][0] * PF_mat[15][10];
	P_mat[77] += matrix->values[11][0] * PF_mat[15][11];
	P_mat[66] += matrix->values[11][1] * PF_mat[16][0];
	P_mat[67] += matrix->values[11][1] * PF_mat[16][1];
	P_mat[68] += matrix->values[11][1] * PF_mat[16][2];
	P_mat[69] += matrix->values[11][1] * PF_mat[16][3];
	P_mat[70] += matrix->values[11][1] * PF_mat[16][4];
	P_mat[71] += matrix->values[11][1] * PF_mat[16][5];
	P_mat[72] += matrix->values[11][1] * PF_mat[16][6];
	P_mat[75] += matrix->values[11][1] * PF_mat[16][9];
	P_mat[76] += matrix->values[11][1] * PF_mat[16][10];
	P_mat[77] += matrix->values[11][1] * PF_mat[16][11];
	P_mat[66] += matrix->values[11][2] * PF_mat[17][0];
	P_mat[67] += matrix->values[11][2] * PF_mat[17][1];
	P_mat[68] += matrix->values[11][2] * PF_mat[17][2];
	P_mat[69] += matrix->values[11][2] * PF_mat[17][3];
	P_mat[70] += matrix->values[11][2] * PF_mat[17][4];
	P_mat[71] += matrix->values[11][2] * PF_mat[17][5];
	P_mat[72] += matrix->values[11][2] * PF_mat[17][6];
	P_mat[75] += matrix->values[11][2] * PF_mat[17][9];
	P_mat[76] += matrix->values[11][2] * PF_mat[17][10];
	P_mat[77] += matrix->values[11][2] * PF_mat[17][11];
#else
	P_mat[0]  += matrix->values[0][0] * PF_mat[3][0];
	P_mat[1]  += matrix->values[1][0] * PF_mat[4][0];
	P_mat[2]  += matrix->values[1][0] * PF_mat[4][1];
	P_mat[3]  += matrix->values[2][0] * PF_mat[5][0];
	P_mat[4]  += matrix->values[2][0] * PF_mat[5][1];
	P_mat[5]  += matrix->values[2][0] * PF_mat[5][2];
	P_mat[6]  += matrix->values[3][0] * PF_mat[0][0];
	P_mat[7]  += matrix->values[3][0] * PF_mat[0][1];
	P_mat[8]  += matrix->values[3][0] * PF_mat[0][2];
	P_mat[9]  += matrix->values[3][0] * PF_mat[0][3];
	P_mat[6]  += matrix->values[3][1] * PF_mat[4][0];
	P_mat[7]  += matrix->values[3][1] * PF_mat[4][1];
	P_mat[8]  += matrix->values[3][1] * PF_mat[4][2];
	P_mat[9]  += matrix->values[3][1] * PF_mat[4][3];
	P_mat[6]  += matrix->values[3][2] * PF_mat[11][0];
	P_mat[7]  += matrix->values[3][2] * PF_mat[11][1];
	P_mat[8]  += matrix->values[3][2] * PF_mat[11][2];
	P_mat[9]  += matrix->values[3][2] * PF_mat[11][3];
	P_mat[6]  += matrix->values[3][3] * PF_mat[12][0];
	P_mat[7]  += matrix->values[3][3] * PF_mat[12][1];
	P_mat[8]  += matrix->values[3][3] * PF_mat[12][2];
	P_mat[9]  += matrix->values[3][3] * PF_mat[12][3];
	P_mat[6]  += matrix->values[3][4] * PF_mat[13][0];
	P_mat[7]  += matrix->values[3][4] * PF_mat[13][1];
	P_mat[8]  += matrix->values[3][4] * PF_mat[13][2];
	P_mat[9]  += matrix->values[3][4] * PF_mat[13][3];
	P_mat[6]  += matrix->values[3][5] * PF_mat[14][0];
	P_mat[7]  += matrix->values[3][5] * PF_mat[14][1];
	P_mat[8]  += matrix->values[3][5] * PF_mat[14][2];
	P_mat[9]  += matrix->values[3][5] * PF_mat[14][3];
	P_mat[6]  += matrix->values[3][6] * PF_mat[15][0];
	P_mat[7]  += matrix->values[3][6] * PF_mat[15][1];
	P_mat[8]  += matrix->values[3][6] * PF_mat[15][2];
	P_mat[9]  += matrix->values[3][6] * PF_mat[15][3];
	P_mat[10] += matrix->values[4][0] * PF_mat[1][0];
	P_mat[11] += matrix->values[4][0] * PF_mat[1][1];
	P_mat[12] += matrix->values[4][0] * PF_mat[1][2];
	P_mat[13] += matrix->values[4][0] * PF_mat[1][3];
	P_mat[14] += matrix->values[4][0] * PF_mat[1][4];
	P_mat[10] += matrix->values[4][1] * PF_mat[3][0];
	P_mat[11] += matrix->values[4][1] * PF_mat[3][1];
	P_mat[12] += matrix->values[4][1] * PF_mat[3][2];
	P_mat[13] += matrix->values[4][1] * PF_mat[3][3];
	P_mat[14] += matrix->values[4][1] * PF_mat[3][4];
	P_mat[10] += matrix->values[4][2] * PF_mat[10][0];
	P_mat[11] += matrix->values[4][2] * PF_mat[10][1];
	P_mat[12] += matrix->values[4][2] * PF_mat[10][2];
	P_mat[13] += matrix->values[4][2] * PF_mat[10][3];
	P_mat[14] += matrix->values[4][2] * PF_mat[10][4];
	P_mat[10] += matrix->values[4][3] * PF_mat[12][0];
	P_mat[11] += matrix->values[4][3] * PF_mat[12][1];
	P_mat[12] += matrix->values[4][3] * PF_mat[12][2];
	P_mat[13] += matrix->values[4][3] * PF_mat[12][3];
	P_mat[14] += matrix->values[4][3] * PF_mat[12][4];
	P_mat[10] += matrix->values[4][4] * PF_mat[13][0];
	P_mat[11] += matrix->values[4][4] * PF_mat[13][1];
	P_mat[12] += matrix->values[4][4] * PF_mat[13][2];
	P_mat[13] += matrix->values[4][4] * PF_mat[13][3];
	P_mat[14] += matrix->values[4][4] * PF_mat[13][4];
	P_mat[10] += matrix->values[4][5] * PF_mat[14][0];
	P_mat[11] += matrix->values[4][5] * PF_mat[14][1];
	P_mat[12] += matrix->values[4][5] * PF_mat[14][2];
	P_mat[13] += matrix->values[4][5] * PF_mat[14][3];
	P_mat[14] += matrix->values[4][5] * PF_mat[14][4];
	P_mat[10] += matrix->values[4][6] * PF_mat[15][0];
	P_mat[11] += matrix->values[4][6] * PF_mat[15][1];
	P_mat[12] += matrix->values[4][6] * PF_mat[15][2];
	P_mat[13] += matrix->values[4][6] * PF_mat[15][3];
	P_mat[14] += matrix->values[4][6] * PF_mat[15][4];
	P_mat[15] += matrix->values[5][0] * PF_mat[10][0];
	P_mat[16] += matrix->values[5][0] * PF_mat[10][1];
	P_mat[17] += matrix->values[5][0] * PF_mat[10][2];
	P_mat[18] += matrix->values[5][0] * PF_mat[10][3];
	P_mat[19] += matrix->values[5][0] * PF_mat[10][4];
	P_mat[20] += matrix->values[5][0] * PF_mat[10][5];
	P_mat[15] += matrix->values[5][1] * PF_mat[11][0];
	P_mat[16] += matrix->values[5][1] * PF_mat[11][1];
	P_mat[17] += matrix->values[5][1] * PF_mat[11][2];
	P_mat[18] += matrix->values[5][1] * PF_mat[11][3];
	P_mat[19] += matrix->values[5][1] * PF_mat[11][4];
	P_mat[20] += matrix->values[5][1] * PF_mat[11][5];
	P_mat[15] += matrix->values[5][2] * PF_mat[13][0];
	P_mat[16] += matrix->values[5][2] * PF_mat[13][1];
	P_mat[17] += matrix->values[5][2] * PF_mat[13][2];
	P_mat[18] += matrix->values[5][2] * PF_mat[13][3];
	P_mat[19] += matrix->values[5][2] * PF_mat[13][4];
	P_mat[20] += matrix->values[5][2] * PF_mat[13][5];
	P_mat[15] += matrix->values[5][3] * PF_mat[14][0];
	P_mat[16] += matrix->values[5][3] * PF_mat[14][1];
	P_mat[17] += matrix->values[5][3] * PF_mat[14][2];
	P_mat[18] += matrix->values[5][3] * PF_mat[14][3];
	P_mat[19] += matrix->values[5][3] * PF_mat[14][4];
	P_mat[20] += matrix->values[5][3] * PF_mat[14][5];
	P_mat[15] += matrix->values[5][4] * PF_mat[15][0];
	P_mat[16] += matrix->values[5][4] * PF_mat[15][1];
	P_mat[17] += matrix->values[5][4] * PF_mat[15][2];
	P_mat[18] += matrix->values[5][4] * PF_mat[15][3];
	P_mat[19] += matrix->values[5][4] * PF_mat[15][4];
	P_mat[20] += matrix->values[5][4] * PF_mat[15][5];
	P_mat[21] += matrix->values[6][0] * PF_mat[9][0];
	P_mat[22] += matrix->values[6][0] * PF_mat[9][1];
	P_mat[23] += matrix->values[6][0] * PF_mat[9][2];
	P_mat[24] += matrix->values[6][0] * PF_mat[9][3];
	P_mat[25] += matrix->values[6][0] * PF_mat[9][4];
	P_mat[26] += matrix->values[6][0] * PF_mat[9][5];
	P_mat[27] += matrix->values[6][0] * PF_mat[9][6];
	P_mat[55] += matrix->values[10][0] * PF_mat[11][0];
	P_mat[56] += matrix->values[10][0] * PF_mat[11][1];
	P_mat[57] += matrix->values[10][0] * PF_mat[11][2];
	P_mat[58] += matrix->values[10][0] * PF_mat[11][3];
	P_mat[59] += matrix->values[10][0] * PF_mat[11][4];
	P_mat[60] += matrix->values[10][0] * PF_mat[11][5];
	P_mat[61] += matrix->values[10][0] * PF_mat[11][6];
	P_mat[62] += matrix->values[10][0] * PF_mat[11][7];
	P_mat[63] += matrix->values[10][0] * PF_mat[11][8];
	P_mat[64] += matrix->values[10][0] * PF_mat[11][9];
	P_mat[65] += matrix->values[10][0] * PF_mat[11][10];
	P_mat[55] += matrix->values[10][1] * PF_mat[16][0];
	P_mat[56] += matrix->values[10][1] * PF_mat[16][1];
	P_mat[57] += matrix->values[10][1] * PF_mat[16][2];
	P_mat[58] += matrix->values[10][1] * PF_mat[16][3];
	P_mat[59] += matrix->values[10][1] * PF_mat[16][4];
	P_mat[60] += matrix->values[10][1] * PF_mat[16][5];
	P_mat[61] += matrix->values[10][1] * PF_mat[16][6];
	P_mat[62] += matrix->values[10][1] * PF_mat[16][7];
	P_mat[63] += matrix->values[10][1] * PF_mat[16][8];
	P_mat[64] += matrix->values[10][1] * PF_mat[16][9];
	P_mat[65] += matrix->values[10][1] * PF_mat[16][10];
	P_mat[55] += matrix->values[10][2] * PF_mat[17][0];
	P_mat[56] += matrix->values[10][2] * PF_mat[17][1];
	P_mat[57] += matrix->values[10][2] * PF_mat[17][2];
	P_mat[58] += matrix->values[10][2] * PF_mat[17][3];
	P_mat[59] += matrix->values[10][2] * PF_mat[17][4];
	P_mat[60] += matrix->values[10][2] * PF_mat[17][5];
	P_mat[61] += matrix->values[10][2] * PF_mat[17][6];
	P_mat[62] += matrix->values[10][2] * PF_mat[17][7];
	P_mat[63] += matrix->values[10][2] * PF_mat[17][8];
	P_mat[64] += matrix->values[10][2] * PF_mat[17][9];
	P_mat[65] += matrix->values[10][2] * PF_mat[17][10];
	P_mat[55] += matrix->values[10][3] * PF_mat[18][0];
	P_mat[56] += matrix->values[10][3] * PF_mat[18][1];
	P_mat[57] += matrix->values[10][3] * PF_mat[18][2];
	P_mat[58] += matrix->values[10][3] * PF_mat[18][3];
	P_mat[59] += matrix->values[10][3] * PF_mat[18][4];
	P_mat[60] += matrix->values[10][3] * PF_mat[18][5];
	P_mat[61] += matrix->values[10][3] * PF_mat[18][6];
	P_mat[62] += matrix->values[10][3] * PF_mat[18][7];
	P_mat[63] += matrix->values[10][3] * PF_mat[18][8];
	P_mat[64] += matrix->values[10][3] * PF_mat[18][9];
	P_mat[65] += matrix->values[10][3] * PF_mat[18][10];
	P_mat[66] += matrix->values[11][0] * PF_mat[10][0];
	P_mat[67] += matrix->values[11][0] * PF_mat[10][1];
	P_mat[68] += matrix->values[11][0] * PF_mat[10][2];
	P_mat[69] += matrix->values[11][0] * PF_mat[10][3];
	P_mat[70] += matrix->values[11][0] * PF_mat[10][4];
	P_mat[71] += matrix->values[11][0] * PF_mat[10][5];
	P_mat[72] += matrix->values[11][0] * PF_mat[10][6];
	P_mat[73] += matrix->values[11][0] * PF_mat[10][7];
	P_mat[74] += matrix->values[11][0] * PF_mat[10][8];
	P_mat[75] += matrix->values[11][0] * PF_mat[10][9];
	P_mat[76] += matrix->values[11][0] * PF_mat[10][10];
	P_mat[77] += matrix->values[11][0] * PF_mat[10][11];
	P_mat[66] += matrix->values[11][1] * PF_mat[16][0];
	P_mat[67] += matrix->values[11][1] * PF_mat[16][1];
	P_mat[68] += matrix->values[11][1] * PF_mat[16][2];
	P_mat[69] += matrix->values[11][1] * PF_mat[16][3];
	P_mat[70] += matrix->values[11][1] * PF_mat[16][4];
	P_mat[71] += matrix->values[11][1] * PF_mat[16][5];
	P_mat[72] += matrix->values[11][1] * PF_mat[16][6];
	P_mat[73] += matrix->values[11][1] * PF_mat[16][7];
	P_mat[74] += matrix->values[11][1] * PF_mat[16][8];
	P_mat[75] += matrix->values[11][1] * PF_mat[16][9];
	P_mat[76] += matrix->values[11][1] * PF_mat[16][10];
	P_mat[77] += matrix->values[11][1] * PF_mat[16][11];
	P_mat[66] += matrix->values[11][2] * PF_mat[17][0];
	P_mat[67] += matrix->values[11][2] * PF_mat[17][1];
	P_mat[68] += matrix->values[11][2] * PF_mat[17][2];
	P_mat[69] += matrix->values[11][2] * PF_mat[17][3];
	P_mat[70] += matrix->values[11][2] * PF_mat[17][4];
	P_mat[71] += matrix->values[11][2] * PF_mat[17][5];
	P_mat[72] += matrix->values[11][2] * PF_mat[17][6];
	P_mat[73] += matrix->values[11][2] * PF_mat[17][7];
	P_mat[74] += matrix->values[11][2] * PF_mat[17][8];
	P_mat[75] += matrix->values[11][2] * PF_mat[17][9];
	P_mat[76] += matrix->values[11][2] * PF_mat[17][10];
	P_mat[77] += matrix->values[11][2] * PF_mat[17][11];
	P_mat[66] += matrix->values[11][3] * PF_mat[18][0];
	P_mat[67] += matrix->values[11][3] * PF_mat[18][1];
	P_mat[68] += matrix->values[11][3] * PF_mat[18][2];
	P_mat[69] += matrix->values[11][3] * PF_mat[18][3];
	P_mat[70] += matrix->values[11][3] * PF_mat[18][4];
	P_mat[71] += matrix->values[11][3] * PF_mat[18][5];
	P_mat[72] += matrix->values[11][3] * PF_mat[18][6];
	P_mat[73] += matrix->values[11][3] * PF_mat[18][7];
	P_mat[74] += matrix->values[11][3] * PF_mat[18][8];
	P_mat[75] += matrix->values[11][3] * PF_mat[18][9];
	P_mat[76] += matrix->values[11][3] * PF_mat[18][10];
	P_mat[77] += matrix->values[11][3] * PF_mat[18][11];
	P_mat[78] += matrix->values[12][0] * PF_mat[16][0];
	P_mat[79] += matrix->values[12][0] * PF_mat[16][1];
	P_mat[80] += matrix->values[12][0] * PF_mat[16][2];
	P_mat[81] += matrix->values[12][0] * PF_mat[16][3];
	P_mat[82] += matrix->values[12][0] * PF_mat[16][4];
	P_mat[83] += matrix->values[12][0] * PF_mat[16][5];
	P_mat[84] += matrix->values[12][0] * PF_mat[16][6];
	P_mat[85] += matrix->values[12][0] * PF_mat[16][7];
	P_mat[86] += matrix->values[12][0] * PF_mat[16][8];
	P_mat[87] += matrix->values[12][0] * PF_mat[16][9];
	P_mat[88] += matrix->values[12][0] * PF_mat[16][10];
	P_mat[89] += matrix->values[12][0] * PF_mat[16][11];
	P_mat[90] += matrix->values[12][0] * PF_mat[16][12];
	P_mat[78] += matrix->values[12][1] * PF_mat[17][0];
	P_mat[79] += matrix->values[12][1] * PF_mat[17][1];
	P_mat[80] += matrix->values[12][1] * PF_mat[17][2];
	P_mat[81] += matrix->values[12][1] * PF_mat[17][3];
	P_mat[82] += matrix->values[12][1] * PF_mat[17][4];
	P_mat[83] += matrix->values[12][1] * PF_mat[17][5];
	P_mat[84] += matrix->values[12][1] * PF_mat[17][6];
	P_mat[85] += matrix->values[12][1] * PF_mat[17][7];
	P_mat[86] += matrix->values[12][1] * PF_mat[17][8];
	P_mat[87] += matrix->values[12][1] * PF_mat[17][9];
	P_mat[88] += matrix->values[12][1] * PF_mat[17][10];
	P_mat[89] += matrix->values[12][1] * PF_mat[17][11];
	P_mat[90] += matrix->values[12][1] * PF_mat[17][12];
	P_mat[78] += matrix->values[12][2] * PF_mat[18][0];
	P_mat[79] += matrix->values[12][2] * PF_mat[18][1];
	P_mat[80] += matrix->values[12][2] * PF_mat[18][2];
	P_mat[81] += matrix->values[12][2] * PF_mat[18][3];
	P_mat[82] += matrix->values[12][2] * PF_mat[18][4];
	P_mat[83] += matrix->values[12][2] * PF_mat[18][5];
	P_mat[84] += matrix->values[12][2] * PF_mat[18][6];
	P_mat[85] += matrix->values[12][2] * PF_mat[18][7];
	P_mat[86] += matrix->values[12][2] * PF_mat[18][8];
	P_mat[87] += matrix->values[12][2] * PF_mat[18][9];
	P_mat[88] += matrix->values[12][2] * PF_mat[18][10];
	P_mat[89] += matrix->values[12][2] * PF_mat[18][11];
	P_mat[90] += matrix->values[12][2] * PF_mat[18][12];
#endif
}
#endif

/* This function is used to propagate the covariance matrix in the Kalman Filter, and
 * implements a generic function that can be invoked to update both the 18 and the 9
 * states Kalman Filter. In this case, the sparse matrix is the transition matrix F
 * P(k) = F*P(k-1)*F'
 */
void update_matrix_sparse(unsigned int size, const sparse_matrix_t *matrix, float *P_mat)
{
	float cov_mat [KNSTATES][KNSTATES];
	float PF_mat  [KNSTATES][KNSTATES];

	int i, j, k;

	int index = 0;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j <= i; j++)
		{
			PF_mat  [i][j] = 0.0;
			PF_mat  [j][i] = 0.0;
			cov_mat [i][j] = P_mat[index];
			cov_mat [j][i] = P_mat[index];
			index++;
		}
	}

	for (j = 0; j < size; j++)
	{
		for (k = 0; k < matrix->num_elem_per_row[j]; k++)
		{
			int   col_index = matrix->col_number [j][k];
			float temp      = matrix->values     [j][k];

			for (i = 0; i < size; i++)
			{
				PF_mat[i][j] += cov_mat[col_index][i] * temp;
			}
		}
	}

	float *p_ptr = &P_mat[0];
	for (i = 0; i < size; i++)
	{
		for (j = 0; j <= i; j++)
		{
			p_ptr[j] *= matrix->diagonal[i] * matrix->diagonal[j];
			p_ptr[j] += matrix->diagonal[i] * PF_mat[i][j];
			p_ptr[j] += matrix->diagonal[j] * PF_mat[j][i];
		}

		for (k = 0; k < matrix->num_elem_per_row[i]; k++)
		{
			int col_index = matrix->col_number [i][k];
			float temp    = matrix->values     [i][k];
			for (j = 0; j <= i; j++)
			{
				p_ptr[j] += temp * PF_mat[col_index][j];
			}
		}
		p_ptr += i+1;
	}
}

/* This function performs the Cholesky decomposition of a matrix */
int choleskyDecomp(float *pdMatrixA, int siSize)
{
	/* Auxiliar indices for the different loops */
	int siI, siJ, siK;

	/* Auxiliar variable to sum the contributions to an element's update */
	float dSum;

	/* Pointer to a row of the A or L matrices  */
	float *pdRowI;
	float *pdRowJ;

	/*----------- Input checking --------------*/
	/* Check dimensions */
	if ( siSize < 0 )
	{
		return 0;
	}
	else if (siSize == 0)
	{
		return 1;
	}

	/* Check input pointers */
	if ( pdMatrixA == 0 )
	{
		return 0;
	}

	/*
	 * Input data are valid.
	 * Computations are performed.
	 */

	/***************************************************************************
	 * The objective is to obtain the Cholesky decomposition:                  *
	 *    A = L*Lt                                                             *
	 * where L is a lower triangular matrix and Lt its transpose.              *
	 *                                                                         *
	 * The resulting equations can be solved iteratively. The first row is:    *
	 *    L(1,1) = sqrt(A(1,1))                                                *
	 * The next rows follow from the expressions:                              *
	 *    L(i,j) = (A(i,j) - Sum(k=1,...,j-1){L(i,k)*L(j,k)}) / L(j,j), i > j  *
	 *    L(i,i) = sqrt(A(i,i) - Sum(k=1,...,i-1){L(i,k)*L(i,k)})              *
	 *                                                                         *
	 * Note that the computed values of L are written over A. This is correct, *
	 * since the elements of A used in the above expressions always correspond *
	 * to indices for which L has not been yet computed.                       *
	 *                                                                         *
	 * Also note that in the implementation indices begin by 0 instead of 1.   *
	 *                                                                         *
	 * The way of moving around A and L in the packed format is to perform     *
	 * the correct updates of the index or pointer q that accesses A or L:     *
	 *    L(i,j) => L(i,j+1)  ~  q => q + 1                                    *
	 *    L(i,j) => L(i+1,j)  ~  q => q + i + 1                                *
	 ***************************************************************************/

	/* The first row begins in element 0 */
	pdRowI = &(pdMatrixA[0]);

	/* Solve iteratively row by row (index i) */
	for (siI = 0; siI < siSize; siI = siI + 1)
	{
		/****************************
		 *  Compute L(i,j), j < i)  *
		 * **************************/
		/* First column j is 0. Point pdRowJ to first element. */
		pdRowJ = &(pdMatrixA[0]);
		/* For each row, solve column by column (index j) */
		for (siJ = 0; siJ < siI; siJ = siJ + 1)
		{
			/* Compute (A(i,j) - Sum(k=1,...,j-1){L(i,k)*L(j,k)}) / L(j,j),
			 * in 3 steps:
			 * 1. s = A(i,j) */
			dSum = pdRowI[siJ];

			/* 2. s = s - Sum(k=1,...,j-1){L(i,k)*L(j,k)}) */
			const float * auxI = pdRowI - 1;
			const float * auxJ = pdRowJ - 1;
			for (siK = siJ; siK; siK--)
			{
				dSum = dSum - auxI[siK] * auxJ[siK];
			}

			/* 3. L(i,j) = s / L(j,j).
			 * [The j-th element in row j is the diagonal element L(j,j)] */
			pdRowI[siJ] = dSum / pdRowJ[siJ];

			/* Jump to next row: L(j,0) => L(j+1,0)  ~  q => q + j + 1 */
			pdRowJ = pdRowJ + siJ + 1;

		}/* End loop over columns of L (index j) */

		/********************
		 *  Compute L(i,i), *
		 * ******************/
		/* L(i,i) = sqrt(A(i,i) - Sum(k=1,...,i-1){L(i,k)*L(i,k)}, again
		 * computed in 3 steps:
		 * 1. s = L(i,i) */
		dSum = pdRowI[siI];

		/* 2. s = s - Sum(k=1,...,j-1){L(i,k)*L(i,k)}) *
		 * Remember that pdRowI[siK] = L(i,k), pdRowJ[siK] = L(j,k) */
		for (siK = 0; siK < siI; siK = siK + 1)
		{
			dSum = dSum - pdRowI[siK] * pdRowI[siK];
		}

		/* 3. L(i,i) = sqrt(s).
		 * Checking the partial sum to obtain the square root */
		if ( dSum > MATEPSILON )
		{
			pdRowI[siI] = sqrt(dSum);
		}
		else
		{
			/* The Cholesky decomposition cannot be computed since A is not
			 * positive definite */
			return 0;
		}

		/* Jump to next row: L(i,0) => L(i+1,0)  ~  q => q + i + 1 */
		pdRowI = pdRowI + siI + 1;

	}/* End loop over columns of L (index i) */

	return 1;
}

void triangularSolve(const float *pdMatrixL, const float *pdInputVector, int siSize, float *pdSolutionVector)
{
	/***************************************************************************
	 * The objective is to solve the triangular system:                        *
	 *    L*x = y                                                              *
	 * where L is a lower triangular matrix.                                   *
	 *                                                                         *
	 * The system can be solved iteratively. The first element is:             *
	 *    x(1) = y(1) / L(1,1)                                                 *
	 * The next rows follow from the expression:                               *
	 *    x(i) = (y(i) - Sum(k=1,...,i-1){L(i,k)*x(k)} ) / L(i,i).             *
	 *                                                                         *
	 * Note that the implementation is such that the vector y and x could be   *
	 * the same. The clue is that the elements of y used in the above          *
	 * expression correspond to indices for which x has not been yet computed. *
	 *                                                                         *
	 * Also note that in the implementation indices begin by 0 instead of 1.   *
	 *                                                                         *
	 * The way of moving around A and L in the packed format is to perform     *
	 * the correct updates of the index or pointer q that accesses A or L:     *
	 *    L(i,j) => L(i,j+1)  ~  q => q + 1                                    *
	 *    L(i,j) => L(i+1,j)  ~  q => q + i + 1                                *
	 ***************************************************************************/

	/* Initialise matrix row pointer to the beginning of first row.*/
	const float *pdRowI = &(pdMatrixL[0]);
	const float *pdRowIN;

	/* Loop over elements of the vector.*/
	const int lastI = ((int)(siSize / 2)) * 2;
	int siI;
	for (siI = 0; siI < lastI; siI += 2)
	{
		pdRowIN = pdRowI + siI + 1;

		/* Obtain x(i) in 3 steps:
		 * 1. s = y(i) */
		float dSum00 = 0.0;
		float dSum01 = 0.0;
		float dSum10 = 0.0;
		float dSum11 = 0.0;

		/* 2. s = s - Sum(k=1,...,i-1){L(i,k)*x(k)} *
		 * Remember that pdRowI[siK] = L(i,k) */
		int siK;
		for (siK = siI; siK; siK -= 2)
		{
			dSum00 += pdRowI[siK-1] * pdSolutionVector[siK-1];
			dSum01 += pdRowI[siK-2] * pdSolutionVector[siK-2];
			dSum10 += pdRowIN[siK-1] * pdSolutionVector[siK-1];
			dSum11 += pdRowIN[siK-2] * pdSolutionVector[siK-2];
		}

		/* 3. x(i) = s / L(i,i) */
		pdSolutionVector[siI] = (pdInputVector[siI] - dSum00 - dSum01) / pdRowI[siI];
		dSum10 += pdRowIN[siI] * pdSolutionVector[siI];
		pdSolutionVector[siI+1] = (pdInputVector[siI+1] - dSum10 - dSum11) / pdRowIN[siI+1];

		/* Jump to next row: L(i,0) => L(i+1,0)  ~  q => q + i + 1 */
		pdRowI = pdRowIN + siI + 2;
	}/* End loop over vector elements */

	for (siI = lastI; siI < siSize; siI++)
	{
		/* Obtain x(i) in 3 steps:
		 * 1. s = y(i) */
		float dSum1 = pdInputVector[siI];
		float dSum2 = 0.0;
		float dSum3 = 0.0;
		float dSum4 = 0.0;

		/* 2. s = s - Sum(k=1,...,i-1){L(i,k)*x(k)} *
		 * Remember that pdRowI[siK] = L(i,k) */
		int lastK = ((int)(siI / 4)) * 4;
		int siK = siI;
		for (; siK > lastK; siK--)
		{
			dSum1 -= pdRowI[siK-1] * pdSolutionVector[siK-1];
		}
		for (; siK; siK -= 4)
		{
			dSum1 -= pdRowI[siK-1] * pdSolutionVector[siK-1];
			dSum2 -= pdRowI[siK-2] * pdSolutionVector[siK-2];
			dSum3 -= pdRowI[siK-3] * pdSolutionVector[siK-3];
			dSum4 -= pdRowI[siK-4] * pdSolutionVector[siK-4];
		}

		/* 3. x(i) = s / L(i,i) */
		pdSolutionVector[siI] = (dSum1 + dSum2 +dSum3 + dSum4) / pdRowI[siI];

		/* Jump to next row: L(i,0) => L(i+1,0)  ~  q => q + i + 1 */
		pdRowI = pdRowI + siI + 1;

	}/* End loop over vector elements */
}

void triangularSolveT(const float *pdMatrixL, const float *pdInputVector, int siSize, float *pdSolutionVector)
{

	/* Auxiliar index for the different loops */
	int siI, siK;

	/* Auxiliar variable to make temporary sums. */
	float dSum;

	/* Pointer to a diagonal element of L, L(i,i) */
	const float *pdLii;

	/* Pointer to an element of matrix L , L(i,k) */
	const float *pdLik;

	/* Index of last element in matrix L */
	int siLastIndex = ((siSize + 1) * siSize) / 2 - 1;

	/***************************************************************************
	 * The objective is to solve the triangular system:                        *
	 *    Lt*x = y                                                             *
	 * where L is a lower triangular matrix.                                   *
	 *                                                                         *
	 * The system can be solved iteratively. The first element is:             *
	 *    x(n) = y(n) / L(n,n)                                                 *
	 * The rest of rows follow from the expression, in decreasing order:       *
	 *    x(i) = (y(i) - Sum(k=i+1,...,n){L(k,i)*x(k)} ) / L(i,i).             *
	 *                                                                         *
	 * Note that the implementation is such that the vector y and x could be   *
	 * the same. The clue is that the elements of y used in the above          *
	 * expression correspond to indices for which x has not been yet computed. *
	 *                                                                         *
	 * Also note that in the implementation firts index is 0 instead of 1.     *
	 *                                                                         *
	 * The way of moving around A and L in the packed format is to perform     *
	 * the correct updates of the index or pointer q that accesses A or L:     *
	 *    L(i,j) => L(i,j+1)  ~  q => q + 1                                    *
	 *    L(i,j) => L(i+1,j)  ~  q => q + i + 1                                *
	 ***************************************************************************/

	/* Initialise matrix pointer to the last element L(n,n). */
	pdLii = &(pdMatrixL[siLastIndex]);

	/* Loop over elements of the vector (decreasing order).*/
	for (siI = siSize - 1; siI >= 0; siI = siI - 1)
	{
		/* Obtain x(i) in 3 steps:
		 * 1. s = y(i) */
		dSum = pdInputVector[siI];

		/* 2. s = s - Sum(k=i+1,...,n){L(k,i)*x(k)} *
		 * At the beginning of the loop, pdLik points to L(i,i) */
		pdLik = pdLii;
		for (siK = siI + 1; siK < siSize; siK = siK + 1)
		{
			/* Update pointer: L(i,k-1) => L(i,k)  ~  q => q + k */
			pdLik = pdLik + siK;
			dSum = dSum - (*pdLik) * pdSolutionVector[siK];
		}

		/* 3. x(i) = s / L(i,i) */
		pdSolutionVector[siI] = dSum / (*pdLii);

		/* Update pointer: L(i,i) => L(i-1,i-1)  ~  q => q + -i - 1 */
		pdLii = pdLii - siI - 1;

	}/* End loop over vector elements */
}

void choleskyInv(float *pdMatrixL, int siSize)
{
	/*VARIABLE DECLARATION.*/

	/* Auxiliar indices for the different loops */
	int siI, siJ, siK;

	/* Pointer to a row of input matrix (lower triangular) */
	float *pdLRow;
	float *pdPki;
	float *pdPkj;
	float *pdPjj;

	/* Auxiliar variable to accumulate sums in the inner loops */
	float dSum;

	/*BEGINING OF CODE.*/
	/***************************************************************************
	 * The input is the lower triangular matrix L, which builds the Cholesky   *
	 * decomposition (L*Lt) of a symmetric positive definite matrix A.         *
	 *                                                                         *
	 * The output is the inverse of A = L*Lt, expressed in packed format.      *
	 *                                                                         *
	 * The steps to compute the inverse are the following:                     *
	 * 1.- Compute the inverse of L, which will be the lower triangular matrix *
	 *     P:    Inv(L) = Pt  ==>  Inv(Lt) = P                                 *
	 * 2.- Compute the product Pt*P, which is the inverse of A:                *
	 *       Inv(A) = Inv(L*Lt) = Inv(Lt)*Inv(L) = Inv(L)t*Inv(L) = P*Pt       *
	 *                                                                         *
	 * Note that matrix P is written over L, and that finally Inv(A) = P*Pt    *
	 * is written over P. The sequence of the computations is such that no     *
	 * inconsistency appears.                                                  *
	 *                                                                         *
	 * Also note that the implementation is not CPU optimal, since the inner   *
	 * loops do not have stride 1. Since this function is not critical, no     *
	 * further effort is required to improve it.                               *
	 ***************************************************************************/

	/* (1) Compute lower triangular matrix P, which is the inverse of L.
	 *     The following recursion is used:
	 *     - Beginning
	 *       P(1,1) = 1 / L(1,1)
	 *     - jth-iteration
	 *       P(j,j) = 1 / L(j,j)
	 *       P(i,j) = - Sum(k=j, ..., i-1){L(i,k)*P(k,j)} / L(i,i),  i=j+1,....,n
	 */

	/* Diagonal elements. Initialise to P(0,0) */
	pdPjj = &(pdMatrixL[0]);

	/* Loop over the columns of  */
	for (siJ = 0; siJ < siSize; siJ = siJ + 1)
	{
		/* Initialising row pointer to P(j,j) */
		pdLRow = pdPjj;

		/* P(j,j) = 1 / L(j,j) */
		(*pdPjj) = 1.0 / (*pdLRow);

		/* Set row pointer to L(j+1,0) */
		pdLRow = pdLRow + 1;

		/* Loop over the rows of P */
		for (siI = siJ+1; siI < siSize; siI = siI + 1)
		{
			/* Initialise column pointer to P(j,j) */
			pdPkj = pdPjj;

			/* Compute - Sum(k=j, ..., i-1){L(i,k)*P(k,j)} */
			dSum = 0.0;
			for (siK = siJ; siK < siI; siK = siK + 1)
			{
				dSum = dSum - pdLRow[siK] * *pdPkj;

				/* Update pointer: P(k,j) => P(k+1,j)  ~  q => q + k + 1 */
				pdPkj = pdPkj + siK + 1;
			}

			/* After the loop, pdPkj points to P(i,j). pdLRow points to L(i,0)
			 * Compute P(i,j) = - Sum(k=j, ..., i-1){L(i,k)*P(k,j)} / L(i,i) */
			*pdPkj = dSum / pdLRow[siI];

			/* Update pointer: L(i,0) => L(i+1,0)  ~  q => q + i + 1 */
			pdLRow = pdLRow + siI + 1;

		}/* End loop over i */

		/* Update pointer: P(j,j) => P(j+1,j+1)  ~  q => q + j + 2 */
		pdPjj = pdPjj + siJ + 2;

	}/* End loop over j */


	/* (2) Compute PtP, which is the required inverse of the original A.
	 *     Recursion:
	 *       i-th iteration
	 *         PtP(i,j) = Sum(k=i,...,n){L(k,i)*L(k,j)} ,  j=0,...,i */

	/* Initialise the row pointer to P(0,0) */
	pdLRow = &(pdMatrixL[0]);
	/* Loop over the rows of PtP */
	for (siI = 0; siI < siSize; siI = siI + 1)
	{
		/* Loop over the columns of PtP */
		for (siJ = 0; siJ <= siI; siJ = siJ + 1)
		{
			/* Initialise pointer pdPki to P(i,i) */
			pdPki = pdLRow + siI;
			/* Initialise pointer pdPki to P(i,j) */
			pdPkj = pdLRow + siJ;

			/* Compute Sum(k=i,...,n){L(k,i)*L(k,j)} */
			dSum = 0.0;
			for (siK = siI; siK < siSize; siK = siK + 1)
			{
				dSum = dSum + *pdPki * *pdPkj;

				/* Update pointer: P(k,j) => P(k+1,j)  ~  q => q + k + 1 */
				pdPkj = pdPkj + siK + 1;
				/* Update pointer: P(k,i) => P(k+1,i)  ~  q => q + k + 1 */
				pdPki = pdPki + siK + 1;
			}

			/* PtP(i,j) = Sum(k=i,...,n){L(k,i)*L(k,j)} */
			pdLRow[siJ] = dSum;

		}/* End loop over j */

		/* Update pointer: L(i,0) => L(i+1,0)  ~  q => q + i + 1 */
		pdLRow = pdLRow + siI + 1;

	}/* End loop over i */
}

