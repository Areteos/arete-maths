package calculus.linearAlgebra;

import linearAlgebra.ThomasTridiagonalMatrixSolver;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

class ThomasTridiagonalMatrixSolverTest {
    /*
        2   1   0   0      1
        1   3   1   0      1
        0   1   3   1      1
        0   0   1   2      1

        2   1   0   0      1
        0  5/2  1   0     1/2
        0   0 13/5  1     4/5
        0   0   0 21/13   9/13

        x4 = 3/7
        x3 = 1/7
        x2 = 1/7
        x1 = 3/7
     */
    double[] a1 = new double[]{1, 1, 1};
    double[] b1 = new double[]{2, 3, 3, 2};
    double[] c1 = new double[]{1, 1, 1};
    double[] d1 = new double[]{1, 1, 1, 1};
    double[] x1 = new double[]{3/7d, 1/7d, 1/7d, 3/7d};

    /*
        2  1.5  0      1
       0.5  2  0.5     2
        0  1.5  2      3

        x = -0.1, 0.8, 0.9
     */
    double[] a2 = new double[]{0.5, 1.5};
    double[] b2 = new double[]{2, 2, 2};
    double[] c2 = new double[]{1.5, 0.5};
    double[] d2 = new double[]{1, 2, 3};
    double[] x2 = new double[]{-0.1, 0.8, 0.9};


    /*
        2   1   0      1
        1   3   1      1
        0   1   2      1

        x = 0.5, 0, 0.5
     */
    double[] a3 = new double[]{1, 1};
    double[] b3 = new double[]{2, 3, 2};
    double[] c3 = new double[]{1, 1};
    double[] d3 = new double[]{1, 1, 1};
    double[] x3 = new double[]{0.5, 0, 0.5};

    /*
        2   2   0      1  // first row fails diagonal dominance
        1   3   1      1
        0   1   2      1
     */
    double[] a4 = new double[]{1, 1};
    double[] b4 = new double[]{2, 3, 2};
    double[] c4 = new double[]{2, 1};
    double[] d4 = new double[]{1, 1, 1};

    /*
        2   1   0      1
        1   2   1      1  // second row fails diagonal dominance
        0   1   2      1
     */
    double[] a5 = new double[]{1, 1};
    double[] b5 = new double[]{2, 2, 2};
    double[] c5 = new double[]{1, 1};
    double[] d5 = new double[]{1, 1, 1};

    /*
        2   1   0      1
        1   3   1      1
        0   1   0      1  // third row fails diagonal dominance
     */
    double[] a6 = new double[]{1, 1};
    double[] b6 = new double[]{2, 3, 0};
    double[] c6 = new double[]{1, 1};
    double[] d6 = new double[]{1, 1, 1};


    @BeforeEach
    void setUp() {
    }

    @AfterEach
    void tearDown() {
    }

    @Test
    void solve() {
        assertArrayEquals(x1, ThomasTridiagonalMatrixSolver.solve(a1, b1, c1, d1, true), 1e-10);
        assertArrayEquals(x2, ThomasTridiagonalMatrixSolver.solve(a2, b2, c2, d2, true), 1e-10);
        assertArrayEquals(x3, ThomasTridiagonalMatrixSolver.solve(a3, b3, c3, d3, true), 1e-10);
        assertThrows(AssertionError.class, () -> ThomasTridiagonalMatrixSolver.solve(a4, b4, c4, d4, true));
        assertThrows(AssertionError.class, () -> ThomasTridiagonalMatrixSolver.solve(a5, b5, c5, d5, true));
        assertThrows(AssertionError.class, () -> ThomasTridiagonalMatrixSolver.solve(a6, b6, c6, d6, true));
    }
}