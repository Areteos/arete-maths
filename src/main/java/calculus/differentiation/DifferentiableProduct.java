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

package calculus.differentiation;

import calculus.differentiation.functionTypes.Constant;
import net.jcip.annotations.Immutable;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static calculus.differentiation.functionTypes.Constant.ZERO;
import static functions.MathsFunctions.roundToSignificantFigures;
import static functions.MiscFunctions.stripTrailingZeros;

@Immutable
public class DifferentiableProduct extends DifferentiableFunction {
    final List<DifferentiableFunction> factors;

    private DifferentiableProduct(double coefficient, List<DifferentiableFunction> factors) {
        super(null, coefficient);
        this.factors = factors;
    }

    /**
     * Factory method for generating DifferentiableProduct DifferentiableFunctions. Note that this method is NOT guaranteed to
     * return a DifferentiableProduct: if the values passed in can be represented more simply within this library, then they will be. The
     * only guarantee on the return type is that it will be some kind of DifferentiableFunction, and that it will be the
     * simplest available representation of the input. Note that every sum in the input will be fully expanded
     * within this method, making it extremely likely that, if there are any sums given as factors, the output will itself
     * be a DifferentiableSum.
     * @param inputFactors The factors of the product
     * @return A DifferentiableFunction representing the product of the input factors
     */
    public static @NotNull DifferentiableFunction construct(DifferentiableFunction... inputFactors) {
        return construct(1, Arrays.asList(inputFactors));
    }

    /**
     * Factory method for generating DifferentiableProduct DifferentiableFunctions. Note that this method is NOT guaranteed to
     * return a DifferentiableProduct: if the values passed in can be represented more simply within this library, then they will be. The
     * only guarantee on the return type is that it will be some kind of DifferentiableFunction, and that it will be the
     * simplest available representation of the input. Note that every sum in the input will be fully expanded
     * within this method, making it extremely likely that, if there are any sums given as factors, the output will itself
     * be a DifferentiableSum.
     * @param coefficient A known coefficient of the entire product
     * @param inputFactors The factors of the product
     * @return A DifferentiableFunction representing the product of the input factors, multiplied by the coefficient
     */
    public static @NotNull DifferentiableFunction construct(double coefficient, DifferentiableFunction... inputFactors) {
        return construct(coefficient, Arrays.asList(inputFactors));
    }

    /**
     * Factory method for generating DifferentiableProduct DifferentiableFunctions. Note that this method is NOT guaranteed to
     * return a DifferentiableProduct: if the values passed in can be represented more simply within this library, then they will be. The
     * only guarantee on the return type is that it will be some kind of DifferentiableFunction, and that it will be the
     * simplest available representation of the input. Note that every sum in the input will be fully expanded
     * within this method, making it extremely likely that, if there are any sums given as factors, the output will itself
     * be a DifferentiableSum.
     * @param coefficient A known coefficient of the entire product
     * @param inputFactors The factors of the product
     * @return A DifferentiableFunction representing the product of the input factors, multiplied by the coefficient
     */
    public static @NotNull DifferentiableFunction construct(double coefficient, List<DifferentiableFunction> inputFactors) {
        List<DifferentiableFunction> factors = new ArrayList<>(inputFactors);

        if (coefficient == 0) {
            return ZERO;
        }
        if (factors.size() == 0) {
            throw new RuntimeException();
        } else if (factors.size() == 1) {
            return factors.get(0);
        }


        Map<Class<? extends DifferentiableFunction>, List<DifferentiableFunction>> aggregatedFactors = new HashMap<>();
        List<DifferentiableSum> sumFactors = new ArrayList<>();

        for (int i = 0; i < factors.size();) {
            DifferentiableFunction factor = factors.get(i++);
            if (ZERO.equals(factor)) {
                return ZERO;  // Zero times anything is zero
            } else {
                coefficient *= factor.getCoefficient();
                factor = factor.withCoefficient(1);

                if (factor instanceof DifferentiableProduct) {
                    factors.addAll(((DifferentiableProduct) factor).factors);
                } else if (factor instanceof DifferentiableSum sumFactor) {
                    sumFactors.add(sumFactor);
                } else if (!(factor instanceof Constant)) {
                    aggregateTerm(aggregatedFactors, factor);
                }
            }
        }

        List<DifferentiableFunction> finalFactors = new ArrayList<>();
        aggregatedFactors.values().forEach(finalFactors::addAll);

        if (!sumFactors.isEmpty()) {
            DifferentiableFunction product = sumFactors.remove(0);
            for (DifferentiableFunction finalFactor : finalFactors) {
                product = product.multiply(finalFactor);
            }
            while (!sumFactors.isEmpty()) {
                product = sumFactors.remove(0).multiply(product);
            }
            return product.multiply(coefficient);
        }

        if (finalFactors.size() == 1) {
            return finalFactors.get(0).withCoefficient(coefficient);
        } if (finalFactors.isEmpty()) {
            return new Constant(coefficient);
        }

        return new DifferentiableProduct(coefficient, finalFactors);
    }

    /**
     * Given a Map keyed by DifferentiableFunction Type to a corresponding List of DifferentiableFunctions of that type,
     * and a factor to add; cleanly and statically add the factor to the Map. This means finding its Type in the map,
     * iterating through the corresponding list until a factor is found that it multiplies with cleanly, removing that
     * factor from the list, and then calling this function again on the same (modified) map to insert the result of the
     * clean multiplication. If no clean multiplications can be performed within the corresponding list, the new factor is
     * simply added to the list without recursion.
     * @param aggregatedFactors A Map containing Lists, each of a particular Type of DifferentiableFunction. The Map is
     *                        keyed by corresponding Types
     * @param factor The factor to be inserted into the map
     * @see DifferentiableProduct#construct(double, List)
     */
    private static void aggregateTerm(Map<Class<? extends DifferentiableFunction>, List<DifferentiableFunction>> aggregatedFactors, DifferentiableFunction factor) {
        List<DifferentiableFunction> aggregatedFactorsOfThisType =
                aggregatedFactors.computeIfAbsent(factor.getClass(), k-> new ArrayList<>());
        for (int i = 0; i < aggregatedFactorsOfThisType.size(); i++) {
            DifferentiableFunction extantTerm = aggregatedFactorsOfThisType.get(i);
            Optional<DifferentiableFunction> product = extantTerm.multiplyCleanly(factor);
            if (product.isPresent()) {
                aggregatedFactorsOfThisType.remove(extantTerm);
                aggregateTerm(aggregatedFactors, product.orElseThrow());
                return;
            }
        }
        aggregatedFactorsOfThisType.add(factor);
    }

    /**
     * The strategy to differentiating a term will be product rule: split off the first factor from the term,
     * and find its derivative. Multiply that derivative with all the other factors, undifferentiated, to get
     * the first term in the total resulting derivative. Then, simply bundle all the other factors into a single
     * term, and find their derivative, then multiply that derivative by the first factor to get the second term.
     * This is a recursive operation.
     *
     * @return The first derivative of the entire term
     */
    @Override
    protected DifferentiableFunction differentiateExplicitly() {

        if (factors.size() == 1) {
            throw new RuntimeException();
        }

        DifferentiableFunction firstFactor = factors.get(0);
        DifferentiableFunction otherFactors = construct(getCoefficient(), factors.subList(1, factors.size()));

        // First term
        DifferentiableFunction firstTerm = construct(1,
                List.of(firstFactor.differentiate(), otherFactors));
        // Second term
        DifferentiableFunction secondTerm = construct(1,
                List.of(firstFactor, otherFactors.differentiate()));

        // Return the sum of the first and second terms
        return DifferentiableSum.construct(getCoefficient(), List.of(firstTerm, secondTerm));
    }

    @Override
    public Optional<DifferentiableFunction> divideCleanly(DifferentiableFunction divisor) {

        if (divisor instanceof DifferentiableProduct productDivisor) {
            DifferentiableFunction currentQuotient = this;
            for (DifferentiableFunction factor : productDivisor.factors) {
                Optional<DifferentiableFunction> result = currentQuotient.divideCleanly(factor);
                if (result.isEmpty()) {
                    return Optional.empty();
                } else {
                    currentQuotient = result.orElseThrow();
                }
            }
            return Optional.of(currentQuotient);


        } else {
            for (DifferentiableFunction factor : factors) {
                Optional<DifferentiableFunction> quotient = factor.divideCleanly(divisor);
                if (quotient.isPresent()) {
                    List<DifferentiableFunction> newFactors = new ArrayList<>(factors);
                    newFactors.remove(factor);
                    newFactors.add(quotient.orElseThrow());
                    return Optional.of(construct(getCoefficient(), newFactors));
                }
            }
        }

        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> multiplyCleanly(DifferentiableFunction multiplicand) {

        if (multiplicand instanceof DifferentiableProduct productMultiplicand) {
            DifferentiableFunction currentProduct = this;
            for (DifferentiableFunction factor : productMultiplicand.factors) {
                Optional<DifferentiableFunction> result = currentProduct.multiplyCleanly(factor);
                if (result.isEmpty()) {
                    return Optional.empty();
                } else {
                    currentProduct = result.orElseThrow();
                }
            }
            return Optional.of(currentProduct);

        } else {
            for (DifferentiableFunction factor : factors) {
                Optional<DifferentiableFunction> product = factor.multiplyCleanly(multiplicand);
                if (product.isPresent()) {
                    List<DifferentiableFunction> newFactors = new ArrayList<>(factors);
                    newFactors.remove(factor);
                    newFactors.add(product.orElseThrow());
                    return Optional.of(construct(getCoefficient(), newFactors));
                }
            }
        }

        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> subtractCleanly(DifferentiableFunction subtrahend) {
        if (Objects.equals(withCoefficient(1), subtrahend.withCoefficient(1))) {
            return Optional.of(withCoefficient(this.getCoefficient() - subtrahend.getCoefficient()));
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> addCleanly(DifferentiableFunction addend) {
        if (Objects.equals(withCoefficient(1), addend.withCoefficient(1))) {
            return Optional.of(withCoefficient(this.getCoefficient() + addend.getCoefficient()));
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> exponentiateCleanly(double exponent) {
        List<DifferentiableFunction> exponentiatedFactors = new ArrayList<>();
        for (DifferentiableFunction factor : factors) {
            Optional<DifferentiableFunction> exponentiatedFactor = factor.exponentiateCleanly(exponent);
            if (exponentiatedFactor.isEmpty()) {
                return Optional.empty();
            } else {
                exponentiatedFactors.add(exponentiatedFactor.orElseThrow());
            }
        }
        return Optional.of(construct(getCoefficient(), exponentiatedFactors));
    }

    @Override
    public @NotNull String toString(int significantFigures) {
        StringBuilder string = new StringBuilder();
        string.append("(%s)".formatted(factors.get(0).toString(significantFigures)));
        for (DifferentiableFunction factor : factors.subList(1, factors.size())) {
            string.append("*(%s)".formatted(factor.toString(significantFigures)));
        }
        return getCoefficient()==1? string.toString() : "%s%s".formatted(stripTrailingZeros(roundToSignificantFigures(getCoefficient(),significantFigures)), string.toString());
    }

    @Override
    public @NotNull DifferentiableFunction withCoefficient(double coefficient) {
        if (coefficient == getCoefficient()) {
            return this;
        }
        return construct(coefficient, factors);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        DifferentiableProduct differentiableProduct = (DifferentiableProduct) o;

        if (Double.compare(differentiableProduct.getCoefficient(), getCoefficient()) != 0) {
            return false;
        }

        List<DifferentiableFunction> otherFactors = new ArrayList<>(differentiableProduct.getFactors());
        if (factors.size() != otherFactors.size()) {
            return false;
        }

        for (DifferentiableFunction factor : factors) {
            if (!otherFactors.remove(factor)) {
                return false;
            }
        }
        return otherFactors.isEmpty();

    }

    @Override
    public int hashCode() {
        List<Integer> hashCodeList = new ArrayList<>();
        for (DifferentiableFunction factor : factors) {
            hashCodeList.add(factor.hashCode());
        }
        Collections.sort(hashCodeList);

        return Objects.hash("*", getCoefficient(), hashCodeList);
    }

    @Override
    public DifferentiableFunction withInputFunction(DifferentiableFunction newInputFunction) {
        return compose(newInputFunction);
    }

    @Override
    public DifferentiableFunction compose(DifferentiableFunction innerFunction) {
        List<DifferentiableFunction> composedFactors = new ArrayList<>(factors.size());
        for (DifferentiableFunction factor : factors) {
            composedFactors.add(factor.compose(innerFunction));
        }

        return construct(getCoefficient(), composedFactors);
    }

    @Override
    public double applyAsDouble(double input) {
        double product = getCoefficient();
        for (DifferentiableFunction factor : factors) {
            product *= factor.applyAsDouble(input);
        }
        return product;
    }

    @Override
    public List<DifferentiableFunction> getFactors() {
        return List.copyOf(factors);
    }

    @Override
    public List<DifferentiableFunction> commonFactors(List<DifferentiableFunction> potentialCommonFactors) {
        List<DifferentiableFunction> establishedCommonFactors = new ArrayList<>();
        for (DifferentiableFunction factor : factors) {
            establishedCommonFactors.addAll(factor.commonFactors(potentialCommonFactors));
        }
        return establishedCommonFactors;
    }
}
