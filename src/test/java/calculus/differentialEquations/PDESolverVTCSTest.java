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

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.function.DoubleBinaryOperator;

import static java.lang.Math.*;
import static org.junit.jupiter.api.Assertions.assertEquals;


class PDESolverVTCSTest {
    static class PDEWithBoundariesAndAnalyticSolution implements DoubleBinaryOperator {
        final DoubleBinaryOperator analyticSolution;
        final DoubleBinaryOperator A;
        final DoubleBinaryOperator B;
        final DoubleBinaryOperator C;

        final BoundaryCondition lowerBound;

        final BoundaryCondition upperBound;

        PDEWithBoundariesAndAnalyticSolution(DoubleBinaryOperator analyticSolution, DoubleBinaryOperator a, DoubleBinaryOperator b, DoubleBinaryOperator c, BoundaryCondition lowerBound, BoundaryCondition upperBound) {
            this.analyticSolution = analyticSolution;
            A = a;
            B = b;
            C = c;
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;
        }


        @Override
        public double applyAsDouble(double t, double x) {
            return analyticSolution.applyAsDouble(t, x);
        }
    }


    final List<PDEWithBoundariesAndAnalyticSolution> problems = new ArrayList<>();

    @BeforeEach
    void setUp() {
        problems.clear();

//         u_t = k*u_xx, boundaries equal 0 at x,L, initial conditions of u = 6*sin(pi*x / L)
        record KL(double k, double L){}
        List<KL> functionConditions = List.of(
                new KL(3, 4),
                new KL(0, 1),
                new KL(10, 10)
        );

        for (KL functionCondition : functionConditions) {
            final double k = functionCondition.k;
            final double L = functionCondition.L;
            problems.add(new PDEWithBoundariesAndAnalyticSolution(
                    (t,x) -> 6 * sin(PI * x / L) * exp(-k * pow(PI / L, 2) * t),
                    (t, x) -> k,
                    (t, x) -> 0,
                    (t, x) -> 0,
                    BoundaryCondition.dirichlet(0, x->0),
                    BoundaryCondition.dirichlet(L, x->0)
            ));
        }

        problems.add(new PDEWithBoundariesAndAnalyticSolution(
                (t, x) -> {
                    double sum = 0;
                    for (double n = 1; n < 100; n++) {
                        sum += 8*n*PI * (1 - pow(-1, n) / sqrt(E)) * exp(t * (-4*n*n*PI*PI - 1)/2) * sin(n * PI * x) / (1 + 4*n*n*PI*PI);
                    }
                    return sum * exp(x/2);
                },
                (t, x) -> 2,
                (t, x) -> -2,
                (t, x) -> 1,
                BoundaryCondition.dirichlet(0, t->0),
                BoundaryCondition.dirichlet(1, t->0)
        ));
    }


    @Test
    void test() {
        int numSpaceSteps = 1000;
        double endTime = 10;
        double timeStepSize = 0.001;

        for (double i : new double[]{0, 0.5, 1}) {

            for (PDEWithBoundariesAndAnalyticSolution problem : problems) {
                double range = problem.upperBound.location() - problem.lowerBound.location();
                double nominalSpaceStep = range / (numSpaceSteps - 0.5);

                PDESolver solver = new PDESolverVTCS(i,
                        x -> problem.analyticSolution.applyAsDouble(0, x),
                        problem.A, problem.B, problem.C,
                        problem.lowerBound, problem.upperBound,
                        nominalSpaceStep, timeStepSize
                );

                double worstError = -1;
                for (double t = 0; t <= endTime; t += timeStepSize) {
    //                double finalT = t;
    //                Map<String, Function<Double, Double>> functionMap = Map.of("Analytic solution", x-> problem.analyticSolution.applyAsDouble(finalT, x), "Approximate solution", x -> solver.applyAsDouble(finalT, x));
    //                functions.Visualiser.drawLineGraph("Analytic vs Approximate solutions at t="+t, functionMap, problem.lowerBoundLocation, problem.upperBoundLocation, numSpaceSteps, false, false, "", "");

                    for (double x = problem.lowerBound.location(); x <= problem.upperBound.location(); x += nominalSpaceStep) {
                        double analyticSolution = problem.analyticSolution.applyAsDouble(t, x);
                        double approximateSolution = solver.applyAsDouble(t, x);
                        double error = Math.abs(analyticSolution - approximateSolution);
                        if (error > worstError) worstError = error;
//                        assertEquals(analyticSolution, approximateSolution, 1e-2, "Unacceptable inaccuracy at t=%s, x=%s with implicit ratio of %s".formatted(t, x, i));
                    }
                }
                System.out.printf("Worst error for implicit ratio of %s was %s%n", i, worstError);
            }
        }
    }
}