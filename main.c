#include <stdio.h>
#include <math.h>
#include "src/matrix.c"
int main(int argc, char **argv)
{
    dComplex exp = {0, 2 * PI / 3};
    dComplex base = {E, 0};
    dComplex result = cPow(base, exp);
    // printc("%6.4f", result, true);
    dComplex dataBlock[] = {
        {1, 0},
        {4, 0},
        {7, 0},
        {2, 0},
        {5, 0},
        {8, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
        {9, 0},
        {2, 0},
        {6, 0},
        {3, 0},
    };
    Matrix matrix1 = mInit(4, 4);
    matrix1.pdata = dataBlock;
    printm(matrix1, "%2.1f ", true);
    printf("\n");
    // singleRowElim(matrix1, 1, 0, 0);
    Matrix matrix2 = gaussElim(matrix1);
    printm(matrix2, "%2.1f ", true);
    Matrix matrix3 = reducedRowElim(matrix2);
    printf("\n");
    printm(matrix3, "%2.1f ", true);
    // printm(matrix3, "%2.1f ", true);
    // pIndex matrix1_index = mapIndex(matrix1);
    // dComplex trace = mTrace(matrix1);
    // rowExchange(matrix1, 1, 0);
    // printm(matrix1, "%1.1f ");
    // printf("\n");
    // Matrix product = mProd(matrix1, matrix1);
    // printm(product, "%1.1f ");
    // printf("\n");
    // mConj(product);
    // printm(product, "%1.1f ");
    // printf("\n");
    // printm(mHermitian(product), "%1.1f ");
    return 0;
}