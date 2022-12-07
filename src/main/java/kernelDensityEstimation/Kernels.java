package kernelDensityEstimation;

import calculus.differentiation.DifferentiableFunction;
import calculus.differentiation.functionTypes.Constant;
import calculus.differentiation.functionTypes.Monomial;
import calculus.differentiation.functionTypes.NaturalExponent;

public abstract class Kernels {
    public static DifferentiableFunction gaussian(double bandwidth) {
        return construct(bandwidth, NaturalExponent.getGaussianDistribution(1, 0));
    }

    public static DifferentiableFunction construct(double bandwidth, DifferentiableFunction function) {
        double inverseBandwidth = 1 / bandwidth;
        function = function.compose(Monomial.construct(inverseBandwidth, 1, null)).multiply(new Constant(inverseBandwidth));
        return function;
    }
}
