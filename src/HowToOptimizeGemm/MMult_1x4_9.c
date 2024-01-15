/* Create macros so that the matrices are stored in column-major order */

#define A(i, j) a[(j) * lda + (i)]
#define B(i, j) b[(j) * ldb + (i)]
#define C(i, j) c[(j) * ldc + (i)]

/* Routine for computing C = A * B + C */

void AddDot(int, double *, int, double *, double *);
void AddDot1x4(int, double *, int, double *, int, double *, int);

void MY_MMult(int m, int n, int k, double *a, int lda,
              double *b, int ldb,
              double *c, int ldc)
{
  int i, j;

  for (j = 0; j < n; j+=4)
  { /* Loop over the rows of C */
    for (i = 0; i < m; i++)
    { /* Loop over the columns of C */
      AddDot1x4(k, &A(i, 0), lda, &B(0, j), ldb, &C(i, j), ldc);
    }
  }
}

void AddDot1x4(int k, double *a, int lda, double *b, int ldb, double *c, int ldc)
{
    int p = 0;
    register double c_00_reg, c_01_reg, c_02_reg, c_03_reg, a_0p_reg;
    double *bp0_pntr, *bp1_pntr, *bp2_pntr, *bp3_pntr;

    c_00_reg = 0.0;
    c_01_reg = 0.0;
    c_02_reg = 0.0;
    c_03_reg = 0.0;

    bp0_pntr = &B(0, 0);
    bp1_pntr = &B(0, 1);
    bp2_pntr = &B(0, 2);
    bp3_pntr = &B(0, 3);

    for (p = 0; p < k; p += 4) {
        a_0p_reg = A(0, p);

        c_00_reg += a_0p_reg * *(bp0_pntr);
        c_01_reg += a_0p_reg * *(bp1_pntr);
        c_02_reg += a_0p_reg * *(bp2_pntr);
        c_03_reg += a_0p_reg * *(bp3_pntr);

        a_0p_reg = A(0, p + 1);

        c_00_reg += a_0p_reg * *(bp0_pntr + 1);
        c_01_reg += a_0p_reg * *(bp1_pntr + 1);
        c_02_reg += a_0p_reg * *(bp2_pntr + 1);
        c_03_reg += a_0p_reg * *(bp3_pntr + 1);

        a_0p_reg = A(0, p + 2);

        c_00_reg += a_0p_reg * *(bp0_pntr + 2);
        c_01_reg += a_0p_reg * *(bp1_pntr + 2);
        c_02_reg += a_0p_reg * *(bp2_pntr + 2);
        c_03_reg += a_0p_reg * *(bp3_pntr + 2);

        a_0p_reg = A(0, p + 3);

        c_00_reg += a_0p_reg * *(bp0_pntr + 3);
        c_01_reg += a_0p_reg * *(bp1_pntr + 3);
        c_02_reg += a_0p_reg * *(bp2_pntr + 3);
        c_03_reg += a_0p_reg * *(bp3_pntr + 3);   

        bp0_pntr += 4;
        bp1_pntr += 4;
        bp2_pntr += 4;
        bp3_pntr += 4;

    }

    C(0, 0) += c_00_reg;
    C(0, 1) += c_01_reg;
    C(0, 2) += c_02_reg;
    C(0, 3) += c_03_reg;
}
