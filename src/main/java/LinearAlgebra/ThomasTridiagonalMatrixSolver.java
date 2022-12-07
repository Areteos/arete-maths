package LinearAlgebra;

import java.util.Arrays;
import java.util.InputMismatchException;
import java.util.List;

/**
 * Solves a given tridiagonal matrix of linear equations using the Thomas method. Each row n represents an equation of the
 * form:
 * <p>
 *     &emsp a<sub>n</sub>x<sub>n-1</sub> + b<sub>n</sub>x<sub>n</sub> + c<sub>n</sub>x<sub>n+1</sub> = d<sub>n</sub>
 * </p>
 * Except the first and last rows, which are of the form:
 * <p>
 *     &emsp b<sub>1</sub>x<sub>1</sub> + c<sub>1</sub>x<sub>2</sub> = d<sub>1</sub>
 * </p>
 * and
 * <p>
 *     &emsp a<sub>N</sub>x<sub>N-1</sub> + b<sub>N</sub>x<sub>N</sub> = d<sub>N</sub>
 * </p>
 * respectively.<p></p>
 * Note that for stability, it is necessary that the absolute value of the middle coefficient of any given row is
 * greater the sum of the absolute values of the other two coefficients, like so:
 * <p>
 *     &emsp |a| + |c| < |b|
 * </p>
 */
public final class ThomasTridiagonalMatrixSolver {
    private ThomasTridiagonalMatrixSolver() {}

    /**
     * Solves a given tridiagonal matrix of linear equations using the Thomas method. Each row n represents an equation of the
     * form:
     * <p>
     *     &emsp a<sub>n</sub>x<sub>n-1</sub> + b<sub>n</sub>x<sub>n</sub> + c<sub>n</sub>x<sub>n+1</sub> = d<sub>n</sub>
     * </p>
     * Except the first and last rows, which are of the form:
     * <p>
     *     &emsp b<sub>1</sub>x<sub>1</sub> + c<sub>1</sub>x<sub>2</sub> = d<sub>1</sub>
     * </p>
     * and
     * <p>
     *     &emsp a<sub>N</sub>x<sub>N-1</sub> + b<sub>N</sub>x<sub>N</sub> = d<sub>N</sub>
     * </p>
     * respectively.<p></p>
     * Note that for stability, it is necessary that the absolute value of the middle coefficient of any given row is
     * greater the sum of the absolute values of the other two coefficients, like so:
     * <p>
     *     &emsp |a| + |c| < |b|
     * </p>
     * If this is checked before calling this method, simply set checkDiagonalDominance to false.
     * @param aValues the coefficients of the x<sub>n-1</sub> term in equation n
     * @param bValues the coefficients of the x<sub>n</sub> term in equation n
     * @param cValues the coefficients of the x<sub>n+1</sub> term in equation n
     * @param dValues the constant right-hand value of equation n
     * @param checkDiagonalDominance Whether to check for the diagonal dominance of the matrix, and so assess whether this method's results can be considered stable
     * @return an array representing the values x<sub>1</sub>, x<sub>2</sub>, ... , x<sub>N</sub>
     */
    public static double[] solve(List<Double> aValues, List<Double> bValues, List<Double> cValues, List<Double> dValues, boolean checkDiagonalDominance) {
        return solve(
                aValues.stream().mapToDouble(d -> d).toArray(),
                bValues.stream().mapToDouble(d -> d).toArray(),
                cValues.stream().mapToDouble(d -> d).toArray(),
                dValues.stream().mapToDouble(d -> d).toArray(),
                checkDiagonalDominance
        );
    }

    /**
     * Solves a given tridiagonal matrix of linear equations using the Thomas method. Each row n represents an equation of the
     * form:
     * <p>
     *     &emsp a<sub>n</sub>x<sub>n-1</sub> + b<sub>n</sub>x<sub>n</sub> + c<sub>n</sub>x<sub>n+1</sub> = d<sub>n</sub>
     * </p>
     * Except the first and last rows, which are of the form:
     * <p>
     *     &emsp b<sub>1</sub>x<sub>1</sub> + c<sub>1</sub>x<sub>2</sub> = d<sub>1</sub>
     * </p>
     * and
     * <p>
     *     &emsp a<sub>N</sub>x<sub>N-1</sub> + b<sub>N</sub>x<sub>N</sub> = d<sub>N</sub>
     * </p>
     * respectively.<p></p>
     * Note that for stability, it is necessary that the absolute value of the middle coefficient of any given row is
     * greater the sum of the absolute values of the other two coefficients, like so:
     * <p>
     *     &emsp |a| + |c| < |b|
     * </p>
     * If this is checked before calling this method, simply set checkDiagonalDominance to false.
     * @param aValues the coefficients of the x<sub>n-1</sub> term in equation n
     * @param bValues the coefficients of the x<sub>n</sub> term in equation n
     * @param cValues the coefficients of the x<sub>n+1</sub> term in equation n
     * @param dValues the constant right-hand value of equation n
     * @param checkDiagonalDominance Whether to check for the diagonal dominance of the matrix, and so assess whether this method's results can be considered stable
     * @return an array representing the values x<sub>1</sub>, x<sub>2</sub>, ... , x<sub>N</sub>
     */
    public static double[] solve(double[] aValues, double[] bValues, double[] cValues, double[] dValues, boolean checkDiagonalDominance) {
        final double[] a = Arrays.copyOf(aValues, aValues.length);
        final double[] b = Arrays.copyOf(bValues, bValues.length);
        final double[] c = Arrays.copyOf(cValues, cValues.length);
        final double[] d = Arrays.copyOf(dValues, dValues.length);

        if (a.length != c.length) {
            throw new InputMismatchException("Off-diagonal coefficients must be equal in number");
        }
        final int n = b.length;

        if (n != d.length) {
            throw new InputMismatchException("Diagonal coefficients and equation constants must be equal in number");
        }
        if (a.length != n - 1) {
            throw new InputMismatchException("Off-diagonal coefficients must be one less in number than diagonal coefficients");
        }

        // Check first and last rows for diagonal dominance
        if (checkDiagonalDominance) {
            assert Math.abs(c[0]) < Math.abs(b[0]) : "row %s failed diagonal dominance condition".formatted(0);
            assert Math.abs(a[n-2]) < Math.abs(b[n-1]) : "row %s failed diagonal dominance condition".formatted(n-1);
        }

        // Forward sweep, gaussian elimination of a-values from top to bottom row
        for (int i = 1; i < n; i++) {

            assert !checkDiagonalDominance || i == n - 1 || Math.abs(a[i - 1]) + Math.abs(c[i]) < Math.abs(b[i]) : "row %s failed diagonal dominance condition".formatted(i);

            final double w = a[i-1] / b[i-1];
            b[i] = b[i] - w * c[i-1];
            d[i] = d[i] - w * d[i-1];
        }

        // Backward sweep, iteratively back-substituting x-values from bottom to top row
        final double[] x = new double[n];
        x[n-1] = d[n-1] / b[n-1];
        for (int i = n-2; i >= 0; i--) {
            x[i] = (d[i] - c[i] * x[i+1]) / b[i];
        }

        return x;
    }
}
