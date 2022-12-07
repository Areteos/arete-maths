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


import calculus.differentiation.DifferentiableFunction;
import calculus.differentiation.DifferentiableSum;
import calculus.differentiation.functionTypes.Constant;
import calculus.differentiation.functionTypes.Monomial;
import types.tuples.Pair;
import functions.IterableFunctions;

import java.util.ArrayList;
import java.util.List;

import static functions.IterableFunctions.unzipped;
import static functions.IterableFunctions.zipped;


public class KernelDensityEstimator {
    public static DifferentiableFunction getKernelDensityEstimator(List<Pair<Double, Double>> samplesAndWeights, DifferentiableFunction kernel) {
        Pair<List<Double>, List<Double>> unzippedSamplesAndWeights = unzipped(samplesAndWeights);
        return getKernelDensityEstimator(unzippedSamplesAndWeights.first(), unzippedSamplesAndWeights.second(), kernel);
    }

    public static DifferentiableFunction getKernelDensityEstimator(Iterable<Double> samples, Iterable<Double> weights, DifferentiableFunction kernel) {

        final double cumulativeWeight = IterableFunctions.getSum(weights);

        List<DifferentiableFunction> kernelList = new ArrayList<>();
        for (Pair<Double, Double> sampleAndWeight : zipped(samples, weights)) {
            double sample = sampleAndWeight.first();
            double weight = sampleAndWeight.second();

            DifferentiableFunction innerFunction = DifferentiableSum.construct(Monomial.X, new Constant(-sample));
            kernelList.add(kernel.compose(innerFunction).multiply(weight));
        }
        return DifferentiableSum.construct(1/cumulativeWeight, kernelList);
    }
}
