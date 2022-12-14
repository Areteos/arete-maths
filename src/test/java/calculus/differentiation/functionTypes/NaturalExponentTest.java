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

package calculus.differentiation.functionTypes;

import functions.MathsFunctions;
import functions.Visualiser;
import org.junit.jupiter.api.Test;
import calculus.differentiation.DifferentiableFunction;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

class NaturalExponentTest {

    @Test
    void construct() {
    }

    @Test
    void differentiateExplicitly() {
    }

    @Test
    void divide() {
    }

    @Test
    void multiply() {
    }

    @Test
    void subtract() {
    }

    @Test
    void add() {
    }

    @Test
    void testToString() {
    }

    @Test
    void changeCoefficient() {
    }

    @Test
    void apply() {
    }

    @Test
    void testEquals() {
    }

    @Test
    void testHashCode() {
    }

    @Test
    void commonFactors() {
    }

    @Test
    void getGaussianDistribution() {
//        DifferentiableFunction currentDerivative = NaturalExponent.getGaussianDistribution(1, 0);
//        List<Visualiser.FunctionInfo> functionInfoList = new ArrayList<>();
//
//        functionInfoList.add(new Visualiser.FunctionInfo("Gaussian", currentDerivative));
//        for (int i = 1; i <= 60; i++) {
//            currentDerivative = currentDerivative.differentiate();
//            DifferentiableFunction factorisedDerivative = currentDerivative.factorise();
//            functionInfoList.add(new Visualiser.FunctionInfo("Derivative %s".formatted(i), factorisedDerivative));
//        }
//
//        Visualiser.drawLineGraph("Gaussian distribution and derivatives", functionInfoList, -5, 5, 10000, false, false, "", "", true);

        for (double i = 0.5; i <= 10; i+=0.5) {
            assertEquals(1, MathsFunctions.integrateApproximately(NaturalExponent.getGaussianDistribution(i, 0), -1000, 1000, 10000), 1e-10);
        }

    }
}