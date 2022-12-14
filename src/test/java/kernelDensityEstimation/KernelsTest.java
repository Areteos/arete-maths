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

package kernelDensityEstimation;

import functions.MathsFunctions;
import org.junit.jupiter.api.Test;
import functions.Visualiser;
import calculus.differentiation.DifferentiableFunction;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

class KernelsTest {

    @Test
    void gaussian() {


        for (double bandwidth : List.of(1e-10, 1e-5, 0.01, 0.1, 0.2, 0.5, 0.75, 1d, 2d, 3d, 4d, 10d)) {
            DifferentiableFunction gaussian = Kernels.gaussian(bandwidth);
            assertEquals(1, MathsFunctions.integrateApproximately(gaussian, -5 * bandwidth, 5*bandwidth, 10000000), 1e-5, "Failed with bandwidth of %s".formatted(bandwidth));
        }
    }

    @Test
    void differentiate() {
        double bandwidth = 1;
        DifferentiableFunction currentDerivative = Kernels.gaussian(bandwidth);
        List<Visualiser.FunctionInfo> functionInfoList = new ArrayList<>();

        functionInfoList.add(new Visualiser.FunctionInfo("Gaussian", currentDerivative));
        for (int i = 1; i <= 3; i++) {
            currentDerivative = currentDerivative.differentiate();
            functionInfoList.add(new Visualiser.FunctionInfo("Derivative %s".formatted(i), currentDerivative));
        }

        Visualiser.drawLineGraph("Gaussian Kernel and derivatives", functionInfoList, -5 * bandwidth, 5 * bandwidth, 10000, false, false, "", "", true);

    }
}