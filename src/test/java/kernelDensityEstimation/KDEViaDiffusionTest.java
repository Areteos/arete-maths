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
import calculus.differentiation.functionTypes.NaturalExponent;

import java.util.List;
import java.util.function.DoubleUnaryOperator;

class KDEViaDiffusionTest {

    @Test
    void improvedSleatherJones_Algorithm1() {
    }

    @Test
    void getDiffusionEstimator() {
        double upperBound = 100;
        double lowerBound = 0;
        double range = upperBound - lowerBound;

//        DoubleUnaryOperator distribution = NaturalExponent.getGaussianDistribution(1, 0.5);
        DoubleUnaryOperator firstPeak = NaturalExponent.getGaussianDistribution(range / 10, (range / 4) + lowerBound);
//        DoubleUnaryOperator firstPeak = NaturalExponent.getGaussianDistribution(range / 10, lowerBound);

        DoubleUnaryOperator secondPeak = NaturalExponent.getGaussianDistribution(range / 10, ((3*range) / 4) + lowerBound);
//        DoubleUnaryOperator secondPeak = NaturalExponent.getGaussianDistribution(range / 10, upperBound);
//
        DoubleUnaryOperator distribution = x -> ((firstPeak.applyAsDouble(x) + secondPeak.applyAsDouble(x)) / 2);

        List<Double> samples = MathsFunctions.generatePoints(distribution, lowerBound, upperBound, 1, 100);



        System.out.println("Done generating samples");

        List<Double> evenWeightings = samples.stream().map(x->1d).toList();

        DoubleUnaryOperator probabilityDensity = KDEViaDiffusion.getDiffusionEstimator(samples, evenWeightings, lowerBound, upperBound);
        DoubleUnaryOperator gaussianDensity = KDEViaDiffusion.getGaussianEstimator(samples, evenWeightings, lowerBound, upperBound);


        Visualiser.drawLineGraph("Diffusion-estimated density",
                List.of(
                        new Visualiser.FunctionInfo("estimated probability density", probabilityDensity),
                        new Visualiser.FunctionInfo("actual probability density", distribution),
                        new Visualiser.FunctionInfo("less-complex estimated density", gaussianDensity)),
                lowerBound, upperBound, 100000, false, false, "Value", "Probability density", true);
//        kernelDensityEstimator.graphProbabilityDensityWithKernelsVisible(1000);

        System.out.println("Integrand of approximation:" + MathsFunctions.integrateApproximately(probabilityDensity, lowerBound, upperBound, 100000));
        System.out.println("Integrand of exact:" + MathsFunctions.integrateApproximately(distribution, lowerBound, upperBound, 100000));
        System.out.println(MathsFunctions.findLocalExtrema(probabilityDensity, lowerBound, upperBound, 100000));

    }
}