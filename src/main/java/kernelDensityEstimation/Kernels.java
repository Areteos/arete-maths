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
