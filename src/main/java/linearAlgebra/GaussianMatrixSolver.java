package linearAlgebra;

import java.util.ArrayList;
import java.util.List;

/**
 * Solves a given enhanced matrix of N linear equations in N unknowns using Gaussian elimination. Each row n represents an equation of the
 * form:
 * <p>
 *     &emsp a<sub>n,1</sub>x<sub>1</sub> + a<sub>n,2</sub>x<sub>2</sub> + ... + a<sub>n,N</sub>x<sub>N</sub> = d<sub>n</sub>
 * </p>
 * Each row n of the given matrix must be of the form:
 * <p>
 *     &emsp a<sub>n,1</sub>  a<sub>n,2</sub>  a<sub>n,3</sub>  ...  a<sub>n,N</sub>  d<sub>n</sub>
 * </p>
 */
public final class GaussianMatrixSolver {
    private GaussianMatrixSolver() {}

    /**
     * Solves a given enhanced matrix of N linear equations in N unknowns using Gaussian elimination. Each row n represents an equation of the
     * form:
     * <p>
     *     &emsp a<sub>n,1</sub>x<sub>1</sub> + a<sub>n,2</sub>x<sub>2</sub> + ... + a<sub>n,N</sub>x<sub>N</sub> = d<sub>n</sub>
     * </p>
     * Each row n of the given matrix must be of the form:
     * <p>
     *     &emsp a<sub>n,1</sub>  a<sub>n,2</sub>  a<sub>n,3</sub>  ...  a<sub>n,N</sub>  d<sub>n</sub>
     * </p>
     * @param matrix The matrix of linear equations to solve
     * @return An array representing the values x<sub>1</sub>, x<sub>2</sub>, ... , x<sub>N</sub>
     */
    public static double[] solve(double[][] matrix) {
        final int n = matrix.length;

        List<Integer> pivotRows = new ArrayList<>();

        for (int col = 0; col < n; col++) {
            double pivotRatio = 0;
            Integer pivotRow = null;
            for (int row = 0; row < n; row++) {
                if (pivotRows.contains(row)) {
                    continue;
                }
                if (col > 0) {  // On all iterations after the first
                    double scaleFactor = matrix[row][col - 1];  // Reduce the current row so that the values below the previous pivot are all 0
                    int previousPivot = pivotRows.get(col-1);
                    for (int i = col-1; i <= n; i++) {  // NOTE that this iteration includes the constant term
                        matrix[row][i] = matrix[row][i] - scaleFactor * matrix[previousPivot][i];
                    }
                }

                double leadingValue = matrix[row][col];
                if (leadingValue == 0) {  // If the leading term in this row is now zero, it cannot be the pivot
                    continue;
                }
                double smallestVal = Double.MAX_VALUE;
                for (int i = 1; i < n; i++) {
                    double absVal = Math.abs(matrix[row][i]);
                    if (absVal != 0 && absVal < smallestVal) {
                        smallestVal = absVal;  // Find the maximum ratio between the leading value and another element in the row
                    }
                }
                double ratio = Math.abs(leadingValue / smallestVal);
                if (ratio > pivotRatio) {  // If this row has the largest ratio so far, it is earmarked as the pivot
                    pivotRow = row;
                    pivotRatio = ratio;
                }
            }
            pivotRows.add(pivotRow);  // Record the eventual pivot row and scale it down by its leading term
            if (pivotRow == null) {
                throw new RuntimeException("Matrix is indeterminate");
            }

            for (int i = n; i >= col; i--) {  // NOTE that this iteration includes the constant term
                matrix[pivotRow][i] = matrix[pivotRow][i] / matrix[pivotRow][col];
            }
        }
        double[] result = new double[n];

        for (int i = n-1; i >= 0; i--) {
            int row = pivotRows.get(i);
            for (int j = n-1; j > i; j--) {
                matrix[row][n] -= matrix[row][j] * result[j];
            }
            result[i] = matrix[row][n];
        }
        return result;
    }
}
