/*
 *    Copyright 2022 Glenn Mamacos
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

package calculus.differentialEquations;

import linearAlgebra.GaussianMatrixSolver;
import linearAlgebra.ThomasTridiagonalMatrixSolver;
import types.tuples.Pair;

import java.util.List;
import java.util.function.BiFunction;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;

/**
 * Numerically solves partial differential equations of the form:<p>
 * &emsp&emsp  f<sub>t</sub> = Af<sub>xx</sub> + Bf<sub>x</sub> + Cf<p>
 * &emsp&emsp where f is an unknown function of t and x, and A, B, and C are known functions of t and x<p>
 * And returns a function that closely approximates the original function f.
 * This is done by using a variably implicit method to produce a tridiagonal matrix of linear
 * equations for each forward time step. The matrix is then solved using the Thomas algorithm. Step sizes in x and t are
 * adjusted to ensure stability of the Thomas algorithm.<p>
 * Values between nodes in x and/or t are calculated by linear interpolation.
 */
public class PDESolverVTCS extends PDESolver {

    private final double implicitRatio;
    private final double complementImplicitRatio;
    private final BiFunction<Double, Double, double[]> abcAtXT;


    /**
     * Numerically solves partial differential equations of the form:<p>
     * &emsp&emsp  f<sub>t</sub> = Af<sub>xx</sub> + Bf<sub>x</sub> + Cf<p>
     * &emsp&emsp where f is an unknown function of t and x, and A, B, and C are known functions of t and x<p>
     * And returns a function that closely approximates the original function f.
     * This is done by using a variably implicit method to produce a tridiagonal matrix of linear
     * equations for each forward time step. The matrix is then solved using the Thomas algorithm. Step sizes in x and t are
     * adjusted to ensure stability of the Thomas algorithm.<p>
     * Values between nodes in x and/or t are calculated by linear interpolation.
     * @param implicitRatio The ratio of implicit to explicit integration to carry out
     * @param functionAtT0 The initial distribution of the solution
     * @param A The time- and space-variable coefficient of the f<sub>xx</sub> term
     * @param B The time- and space-variable coefficient of the f<sub>x</sub> term
     * @param C The time- and space-variable coefficient of the f term
     * @param lowerBound All relevant information about the lower bound of the system
     * @param upperBound All relevant information about the upper bound of the system
     * @param maxSpaceStep The maximum allowable step size in the space domain
     * @param minTimeStep The minimum allowable step size in the time domain
     */
    public PDESolverVTCS(double implicitRatio, DoubleUnaryOperator functionAtT0, DoubleBinaryOperator A, DoubleBinaryOperator B, DoubleBinaryOperator C, BoundaryCondition lowerBound, BoundaryCondition upperBound, double maxSpaceStep, double minTimeStep) {
        super(functionAtT0, lowerBound, upperBound, maxSpaceStep, minTimeStep);

        assert implicitRatio <= 1 && implicitRatio >= 0 : "Implicit ratio must be between 0 and 1 inclusive";

        this.implicitRatio = implicitRatio;
        complementImplicitRatio = 1 - implicitRatio;

        abcAtXT = (t, x) -> {
            double inverseSpaceStep = 1/spaceStep;

            double valueOfA = A.applyAsDouble(t, x);

            double aOverH = valueOfA * inverseSpaceStep;
            double bOver2 = B.applyAsDouble(t, x) / 2;


            return new double[]{
                    (bOver2 - aOverH) * inverseSpaceStep,
                    2 * aOverH * inverseSpaceStep - C.applyAsDouble(t, x),
                    (-aOverH - bOver2) * inverseSpaceStep
            };
        };
    }

    /**
     * Numerically solves partial differential equations of the form:<p>
     * &emsp&emsp  f<sub>t</sub> = Af<sub>xx</sub> + Bf<sub>x</sub> + Cf<p>
     * &emsp&emsp where f is an unknown function of t and x, and A, B, and C are known functions of t and x<p>
     * And returns a function that closely approximates the original function f.
     * This is done by using a variably implicit method to produce a tridiagonal matrix of linear
     * equations for each forward time step. The matrix is then solved using the Thomas algorithm. Step sizes in x and t are
     * adjusted to ensure stability of the Thomas algorithm.<p>
     * Values between nodes in x and/or t are calculated by linear interpolation.
     * @param implicitRatio The ratio of implicit to explicit integration to carry out
     * @param samplesAndWeights A list of point samples and their corresponding weights for an initial distribution
     * @param A The time- and space-variable coefficient of the f<sub>xx</sub> term
     * @param B The time- and space-variable coefficient of the f<sub>x</sub> term
     * @param C The time- and space-variable coefficient of the f term
     * @param lowerBound All relevant information about the lower bound of the system
     * @param upperBound All relevant information about the upper bound of the system
     * @param maxSpaceStep The maximum allowable step size in the space domain
     * @param minTimeStep The minimum allowable step size in the time domain
     */
    public PDESolverVTCS(double implicitRatio, List<Pair<Double, Double>> samplesAndWeights, DoubleBinaryOperator A, DoubleBinaryOperator B, DoubleBinaryOperator C, BoundaryCondition lowerBound, BoundaryCondition upperBound, double maxSpaceStep, double minTimeStep) {
        super(samplesAndWeights, lowerBound, upperBound, maxSpaceStep, minTimeStep);

        assert implicitRatio <= 1 && implicitRatio >= 0 : "Implicit ratio must be between 0 and 1 inclusive";

        this.implicitRatio = implicitRatio;
        complementImplicitRatio = 1 - implicitRatio;

        abcAtXT = (t, x) -> {
            double inverseSpaceStep = 1 / spaceStep;

            double valueOfA = A.applyAsDouble(t, x);

            double aOverH = valueOfA * inverseSpaceStep;
            double bOver2 = B.applyAsDouble(t, x) / 2;


            return new double[]{
                    (bOver2 - aOverH) * inverseSpaceStep,
                    2 * aOverH * inverseSpaceStep - C.applyAsDouble(t, x),
                    (-aOverH - bOver2) * inverseSpaceStep
            };
        };
    }

    @Override
    public void solveUpToTime(double time) {
        while (currentTime <= time) {
            currentTime += timeStep;
            double operatingTime = currentTime + timeStep * implicitRatio;

            int initialOffset = lowerBound.boundaryType() == BoundaryCondition.BoundaryType.DIRICHLET ? 1 : 0;
            int finalOffset = upperBound.boundaryType() == BoundaryCondition.BoundaryType.DIRICHLET ? 1 : 0;

            double[]
                    aArray = new double[nodeLocations.length - initialOffset - finalOffset - 1],
                    bArray = new double[nodeLocations.length - initialOffset - finalOffset],
                    cArray = new double[nodeLocations.length - initialOffset - finalOffset - 1],
                    dArray = new double[nodeLocations.length - initialOffset - finalOffset];
            double[] previousSolutions = nodeSolutionsByTimeStep.get(nodeSolutionsByTimeStep.size()-1);

            double[] currentTimeSolutions = new double[nodeLocations.length];

            // FIRST ROW
            switch (lowerBound.boundaryType()) {
                case DIRICHLET -> {
                    double lowerBoundValue = lowerBound.value(operatingTime);
                    currentTimeSolutions[0] = lowerBoundValue;

                    double nodeLocation = nodeLocations[1];

                    double[] abc = abcAtXT.apply(operatingTime, nodeLocation);
                    double a = abc[0], b = abc[1], c = abc[2];
                    double d =
                            - complementImplicitRatio * a * previousSolutions[0]
                            + (1/timeStep - complementImplicitRatio * b) * previousSolutions[1]
                            - complementImplicitRatio * c * previousSolutions[2]
                            - implicitRatio * a * lowerBoundValue;

                    bArray[0] = 1/timeStep + implicitRatio * b;
                    cArray[0] = implicitRatio * c;
                    dArray[0] = d;
                }
                case NEUMANN -> {
                    double nodeLocation = nodeLocations[0];

                    double[] abc = abcAtXT.apply(operatingTime, nodeLocation);
                    double a = abc[0], b = abc[1], c = abc[2];

                    final double aPlusC = a + c;

                    double d =
                            (1/timeStep - complementImplicitRatio * b) * previousSolutions[0]
                            - complementImplicitRatio * aPlusC * previousSolutions[1]
                            + 2 * a * spaceStep * lowerBound.value(operatingTime);

                    bArray[0] = 1/timeStep + implicitRatio * b;
                    cArray[0] = implicitRatio * aPlusC;
                    dArray[0] = d;
                }
                case ROBIN -> {
                    double nodeLocation = nodeLocations[0];

                    double[] abc = abcAtXT.apply(operatingTime, nodeLocation);
                    double a = abc[0], b = abc[1], c = abc[2];

                    final double fraction1 = 2 * a * spaceStep / lowerBound.order1Coefficient(operatingTime);
                    final double fraction2 = fraction1 * lowerBound.order0Coefficient(operatingTime);

                    final double aPlusC = a + c;

                    double d =
                            (1/timeStep - complementImplicitRatio * (b + fraction2)) * previousSolutions[0]
                            - complementImplicitRatio * aPlusC * previousSolutions[1]
                            + fraction1 * lowerBound.value(operatingTime);

                    bArray[0] = 1/timeStep + implicitRatio * (b + fraction2);
                    cArray[0] = implicitRatio * aPlusC;
                    dArray[0] = d;
                }
            }
            // MIDDLE ROWS
            for (int i=1+initialOffset; i<nodeLocations.length-1-finalOffset; i++) {
                double nodeLocation = nodeLocations[i];

                double[] abc = abcAtXT.apply(currentTime + operatingTime, nodeLocation);
                double a = abc[0], b = abc[1], c = abc[2];

                double d =
                        - complementImplicitRatio * a * previousSolutions[i-1]
                        + (1/timeStep - complementImplicitRatio * b) * previousSolutions[i]
                        - complementImplicitRatio * c * previousSolutions[i + 1];

                aArray[i - initialOffset - 1]   = implicitRatio * a;
                bArray[i - initialOffset]       = 1/timeStep + implicitRatio * b;
                cArray[i - initialOffset]       = implicitRatio * c;
                dArray[i - initialOffset]       = d;
            }
            // LAST ROW
            switch (upperBound.boundaryType()) {
                case DIRICHLET -> {
                    int i = nodeLocations.length - 2;

                    double upperBoundValue = upperBound.value().applyAsDouble(operatingTime);
                    currentTimeSolutions[currentTimeSolutions.length-1] = upperBoundValue;

                    double nodeLocation = nodeLocations[i];

                    double[] abc = abcAtXT.apply(operatingTime, nodeLocation);
                    double a = abc[0], b = abc[1], c = abc[2];
//                    double d = - a * previousSolutions[i-1] + (2/timeStep - b) * previousSolutions[i] - c * previousSolutions[i + 1] - c * upperBoundValue;

                    double d =
                            - complementImplicitRatio * a * previousSolutions[i-1]
                            + (1/timeStep - complementImplicitRatio * b) * previousSolutions[i]
                            - complementImplicitRatio * c * previousSolutions[i+1]
                            - implicitRatio * c * upperBoundValue;

                    aArray[i - initialOffset - 1]   = implicitRatio * a;
                    bArray[i - initialOffset]       = 1/timeStep + implicitRatio * b;
                    dArray[i - initialOffset]       = d;
                }
                case NEUMANN -> {
                    int i = nodeLocations.length - 1;

                    double nodeLocation = nodeLocations[i];

                    double[] abc = abcAtXT.apply(operatingTime, nodeLocation);
                    double a = abc[0], b = abc[1], c = abc[2];

                    final double aPlusC = a + c;

//                    double d = - aPlusC * previousSolutions[i-1] + (2/timeStep - b) * previousSolutions[i] - 4 * c * spaceStep * upperBound.netValue().applyAsDouble(operatingTime);

                    double d =
                            - complementImplicitRatio * aPlusC * previousSolutions[i-1]
                            + (1/timeStep - complementImplicitRatio * b) * previousSolutions[i]
                            - 2 * a * spaceStep * upperBound.value(operatingTime);

                    aArray[i - initialOffset - 1]   = implicitRatio * aPlusC;
                    bArray[i - initialOffset]       = 1/timeStep + implicitRatio * b;
                    dArray[i - initialOffset]       = d;
                }
                case ROBIN -> {
                    int i = nodeLocations.length - 1;


                    double nodeLocation = nodeLocations[0];

                    double[] abc = abcAtXT.apply(operatingTime, nodeLocation);
                    double a = abc[0], b = abc[1], c = abc[2];

                    final double order0Coefficient = upperBound.order0Coefficient().applyAsDouble(operatingTime);
                    final double order1Coefficient = upperBound.order1Coefficient().applyAsDouble(operatingTime);

                    final double fraction1 = 2 * c * spaceStep / order1Coefficient;
                    final double fraction2 = fraction1 * order0Coefficient;

                    final double aPlusC = a + c;

//                    double d = - aPlusC * previousSolutions[i-1] + (2/timeStep - b + fraction2) * previousSolutions[i] - 2 * fraction1 * upperBound.netValue().applyAsDouble(operatingTime);

                    double d =
                            - complementImplicitRatio * aPlusC * previousSolutions[i-1]
                            + (1/timeStep - complementImplicitRatio * (b - fraction2)) * previousSolutions[i]
                            - fraction1 * lowerBound.value(operatingTime);

                    aArray[i - initialOffset - 1]   = implicitRatio * aPlusC;
                    bArray[i - initialOffset]       = 1/timeStep + implicitRatio * (b - fraction2);
                    dArray[i - initialOffset]       = d;
                }
            }
            double[] solution;
            try {
                solution = ThomasTridiagonalMatrixSolver.solve(aArray, bArray, cArray, dArray, true);
            } catch (AssertionError e) {
//                e.printStackTrace();
                double[][] matrix = toMatrix(aArray, bArray, cArray, dArray);
                solution = GaussianMatrixSolver.solve(matrix);
            }

            System.arraycopy(solution, 0, currentTimeSolutions, initialOffset, currentTimeSolutions.length - initialOffset - finalOffset);
            nodeSolutionsByTimeStep.add(currentTimeSolutions);
            timeSteps.add(currentTime);
        }
    }

    public static double[][] toMatrix(double[] a, double[] b, double[] c, double[] d) {
        double[][] matrix = new double[b.length][b.length+1];
        for (int i = 0; i < b.length; i++) {
            if (i > 0) {
                matrix[i][i - 1] = a[i-1];
            }
            matrix[i][i] = b[i];
            if (i < b.length-1) {
                matrix[i][i + 1] = c[i];
            }
            matrix[i][b.length] = d[i];
        }
        return matrix;
    }

}
