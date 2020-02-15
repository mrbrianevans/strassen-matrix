/**
 * This program handles matrices, and has methods for operations to perform on matrices <br>
 * It also includes the <em>Strassen Method of Matrix Multiplication</em>
 * @author Brian Evans
 * @since 2020-02-15
 */
public class Matrix {


    public double[][] array;  // array with elements
    public int n;		       // number of columns, rows



    /**
     * Initialize matrix from the array
     * @param array_ - array
     */
    public Matrix(double[][] array_) {
        array = array_;
        n = array.length;
    }


    /**
     * Initialize the matrix where the number of columns and rows is n
     * @param n_  - number of columns and rows
     */
    public Matrix(int n_) {
        array = new double[n_][n_];
        n = n_;
    }


    /**
     * Perform regular matrix multiplication
     * @param b  - matrix to multiply
     * @return   matrix which is a product
     */
    public Matrix multiply(Matrix b) {
        Matrix product = new Matrix(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                product.array[i][j] = 0;
                for (int k = 0; k < n; k++) {
                    product.array[i][j] += array[i][k] * b.array[k][j];
                }
            }
        }
        return product;
    }


    /**
     * Function to compare two matrices (matrix 'this' and 'b'). Two matrices are
     * said to be equal if they
     * 1) have same size
     * 2) absolute value of the difference between corresponding elements is not
     * more than precision (you can set precision to 0.0000001).
     * @param b - matrix to use for comparison
     * @return true if matrix are equivalent, false otherwise
     */
    public boolean equals(Matrix b) {
        if (this.n != b.n)
            return false;

        double precision = 0.0001;
        for (int i = 0; i < this.n; i++) {
            for (int j = 0; j < this.n; j++) {
                if (Math.abs(this.array[i][j] - b.array[i][j]) > precision)
                    return false;
            }
        }
        return true;
    }

    /**
     * Add two matrices
     * @param b  - matrix to be added
     * @return   - sum
     */
    public Matrix add(Matrix b) {
        double [][] sum = new double[this.n][this.n];

        for (int i = 0; i < this.n; i++) {
            for (int j = 0; j < this.n; j++) {
                sum[i][j] = this.array[i][j] + b.array[i][j];
            }
        }
        return new Matrix(sum);
    }

    /**
     * Subtract one matrix from another
     * @param b  - matrix to be subtracted
     * @return   - difference
     */
    public Matrix minus(Matrix b) {
        double [][] difference = new double[this.n][this.n];

        for (int i = 0; i < this.n; i++) {
            for (int j = 0; j < this.n; j++) {
                difference[i][j] = this.array[i][j] - b.array[i][j];
            }
        }
        return new Matrix(difference);
    }

    /**
     * Return string representing the matrix.
     */
    public String toString() {
        StringBuilder output = new StringBuilder("[");

        output.append("[");
        for (int j = 0; j < this.n-1; j++) {
            output.append(Math.round(this.array[0][j])).append("  ");
        }
        output.append(Math.round(this.array[0][n - 1])).append("]\n");

        for (int i = 1; i < this.n-1; i++) {
            output.append(" [");
            for (int j = 0; j < this.n-1; j++) {
                output.append(Math.round(this.array[i][j])).append("  ");
            }
            output.append(Math.round(this.array[i][n - 1])).append("]\n");
        }
        output.append(" [");
        for (int j = 0; j < this.n-1; j++) {
            output.append(Math.round(this.array[n - 1][j])).append("  ");
        }
        output.append(Math.round(this.array[n - 1][n - 1])).append("]");
        output.append("]");
        return output.toString();
    }



    /**
     * Perform matrix multiplication using Strassen algorithm
     * @param b  - matrix to multiply
     * @return   matrix which is a product
     */
    public Matrix multiplyStrassen(Matrix b) {
        int N = b.n;
        double [][] A = this.array;
        double [][] B = b.array;
        if (N == 2) { // if its a 2x2 matrix, base case
            double[] s = new double[10];
            s[0] = B[0][1] - B[1][1];
            s[1] = A[0][0] + A[0][1];
            s[2] = A[1][0] + A[1][1];
            s[3] = B[1][0] - B[0][0];
            s[4] = A[0][0] + A[1][1];
            s[5] = B[0][0] + B[1][1];
            s[6] = A[0][1] - A[1][1];
            s[7] = B[1][0] + B[1][1];
            s[8] = A[0][0] - A[1][0];
            s[9] = B[0][0] + B[0][1];

            double [] p = new double[7];
            p[0] = A[0][0] * s[0];
            p[1] = s[1] * B[1][1];
            p[2] = s[2] * B[0][0];
            p[3] = s[3] * A[1][1];
            p[4] = s[4] * s[5];
            p[5] = s[6] * s[7];
            p[6] = s[8] * s[9];

            double [][] C = new double[2][2];
            C[0][0] = p[4] + p[3] - p[1] + p[5];
            C[0][1] = p[0] + p[1];
            C[1][0] = p[2] + p[3];
            C[1][1] = p[4] + p[0] - p[2] - p[6];

            return new Matrix(C);
        }
        // if its larger than 2x2, recursively solve

        //this part just splits the two matrices up in quarters
        double [][][][] enormA = new double [2][2][N/2][N/2];
        double [][][][] enormB = new double [2][2][N/2][N/2];
        for (int i = 0; i < N / 2; i++) {
            for (int j = 0; j < N / 2; j++) {
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                        enormA[k][l][i][j] = A[i+(k*N/2)][j+(l*N/2)];
                        enormB[k][l][i][j] = B[i+(k*N/2)][j+(l*N/2)];
                    }
                }
            }
        }

        // convert the double arrays into matrices
        Matrix [][] quadrantsA = new Matrix[2][2];
        Matrix [][] quadrantsB = new Matrix[2][2];

        for (int k = 0; k < 2; k++) {
            for (int l = 0; l < 2; l++) {
                quadrantsA[k][l] = new Matrix(enormA[k][l]);
                quadrantsB[k][l] = new Matrix(enormB[k][l]);
            }
        }

        //create the S matrices
        Matrix [] s = new Matrix[10];
        s[0] = quadrantsB[0][1].minus(quadrantsB[1][1]);
        s[1] = quadrantsA[0][0].add(quadrantsA[0][1]);
        s[2] = quadrantsA[1][0].add(quadrantsA[1][1]);
        s[3] = quadrantsB[1][0].minus(quadrantsB[0][0]);
        s[4] = quadrantsA[0][0].add(quadrantsA[1][1]);
        s[5] = quadrantsB[0][0].add(quadrantsB[1][1]);
        s[6] = quadrantsA[0][1].minus(quadrantsA[1][1]);
        s[7] = quadrantsB[1][0].add(quadrantsB[1][1]);
        s[8] = quadrantsA[0][0].minus(quadrantsA[1][0]);
        s[9] = quadrantsB[0][0].add(quadrantsB[0][1]);

        Matrix [] p = new Matrix[7];
        p[0] = quadrantsA[0][0].multiplyStrassen(s[0]);
        p[1] = s[1].multiplyStrassen(quadrantsB[1][1]);
        p[2] = s[2].multiplyStrassen(quadrantsB[0][0]);
        p[3] = quadrantsA[1][1].multiplyStrassen(s[3]);
        p[4] = s[4].multiplyStrassen(s[5]);
        p[5] = s[6].multiplyStrassen(s[7]);
        p[6] = s[8].multiplyStrassen(s[9]);

        Matrix [][] C = new Matrix[2][2];
        C[0][0] = p[4].add(p[3]).minus(p[1]).add(p[5]);
        C[0][1] = p[0].add(p[1]);
        C[1][0] = p[2].add(p[3]);
        C[1][1] = p[4].add(p[0]).minus(p[2]).minus(p[6]);

        Matrix product = new Matrix(N);

        for (int i = 0; i < N / 2; i++) {
            for (int j = 0; j < N / 2; j++) {
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                        product.array[i+(k*N/2)][j+(l*N/2)] = C[k][l].array[i][j];
                    }
                }
            }
        }

        return product;
    }


}


