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
