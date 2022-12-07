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
import types.tuples.Triple;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import functions.Visualiser;
import calculus.differentiation.functionTypes.NaturalExponent;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.DoubleUnaryOperator;

import static org.junit.jupiter.api.Assertions.assertEquals;

class KernelDensityEstimatorTest {

//    List<Double> samples = new ArrayList<>(Arrays.asList(10d,11d,9d,23d,21d,11d,45d,20d,11d,12d));
//    List<Double> samples = new ArrayList<>(Arrays.asList(-2.1,-1.3, -0.4, 1.9, 5.1, 6.2));

    @BeforeEach
    public void setup() {

    }

    @Test
    void bandwidth() {
    }

    @Test
    void getKernelDensityEstimator() {
        DoubleUnaryOperator trueGauss = NaturalExponent.getGaussianDistribution(10, 0);
        DoubleUnaryOperator bimodalGauss = x -> trueGauss.applyAsDouble(x) + trueGauss.applyAsDouble(x - 5);
        class MultimodalGauss implements DoubleUnaryOperator {
            final List<Double> means;

            public MultimodalGauss(List<Double> means) {
                this.means = means;
            }

            @Override
            public double applyAsDouble(double aDouble) {
                double sum = 0;
                for (Double mean : means) {
                    sum += trueGauss.applyAsDouble(aDouble - mean);
                }
                return sum / means.size();
            }
        }

        double minimum = -100;
        double maximum = 100;
        double range = maximum - minimum;
        int steps = 10000;
        double stepSize = range / (steps+1);

        List<Triple<String, DoubleUnaryOperator, List<Double>>> probabilityDensitiesAndSamples = new ArrayList<>();
        BiConsumer<String, DoubleUnaryOperator> addDistribution = (name, x) -> {
            List<Double> samples = MathsFunctions.generatePoints(x, minimum, maximum, MathsFunctions.findIntervalMinimumAndMaximum(x, minimum, maximum, 10000).second(), 1000);
            double areaUnderTheCurve = MathsFunctions.integrateApproximately(x, minimum, maximum, 10000);
            probabilityDensitiesAndSamples.add(new Triple<>(name, x.andThen(y->y/areaUnderTheCurve), samples));
        };

        addDistribution.accept("Gaussian", trueGauss);
        addDistribution.accept("Bimodal gauss", bimodalGauss);
        addDistribution.accept("sin(x)+1", x -> Math.sin(x) + 1);
        addDistribution.accept("5-x^2", x -> 5 - x*x);
        addDistribution.accept("x+10", x -> x + 10);
        addDistribution.accept("Decimodal gauss", new MultimodalGauss(List.of(-9d, -7d, -5d, -3d, -1d, 1d, 3d, 5d, 7d, 9d)));
        addDistribution.accept("Trimodal gauss", new MultimodalGauss(List.of(-7d, 3d, 8d)));

        List<Double> evenWeightings = probabilityDensitiesAndSamples.get(0).third().stream().map(x->1d).toList();

        for (Triple<String, DoubleUnaryOperator, List<Double>> densityAndSamples : probabilityDensitiesAndSamples) {
            String name = densityAndSamples.first();
            DoubleUnaryOperator density = densityAndSamples.second();
            List<Double> samples = densityAndSamples.third();


            DoubleUnaryOperator kernelDensityEstimator = KDEViaDiffusion.getGaussianEstimator(samples, evenWeightings, minimum, maximum);

            List<Visualiser.FunctionInfo> functionInfoList = List.of(
                    new Visualiser.FunctionInfo("Estimated Density", kernelDensityEstimator),
                    new Visualiser.FunctionInfo("Actual Density", densityAndSamples.second())
            );

            Visualiser.drawLineGraph("Estimated Density for "+densityAndSamples.first(), functionInfoList, minimum, maximum, 10000, false, false, "", "", false);
            assertEquals(1, MathsFunctions.integrateApproximately(kernelDensityEstimator, minimum - range*10, maximum + range*10, 1000000), 1e-10, "Incorrect integral for initial distribution " + densityAndSamples.first());
            for (double x = minimum; x < maximum; x+=stepSize) {
                assertEquals(densityAndSamples.second().applyAsDouble(x), kernelDensityEstimator.applyAsDouble(x), 1e-3, "Error at x=%s for distribution %s".formatted(x, densityAndSamples.first()));
            }
        }
    }
}