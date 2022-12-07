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

import net.jcip.annotations.Immutable;
import org.jetbrains.annotations.NotNull;

import java.util.*;

import static calculus.differentiation.functionTypes.Constant.ZERO;
import static functions.MathsFunctions.roundToSignificantFigures;
import static functions.MiscFunctions.stripTrailingZeros;


/**
 * <p>An immutable class representing a differentiable function of the form:</p>
 * &emsp&emsp f(x) = g<sub>0</sub>(x) + g<sub>1</sub>(x) + ... + g<sub>n</sub>(x)
 */
@Immutable
public class DifferentiableSum extends DifferentiableFunction {
    private final List<DifferentiableFunction> terms;

    private DifferentiableSum(List<DifferentiableFunction> terms) {
        super(null, 1);
        this.terms = terms;  // No need to defensively copy here, since a defensive copy was already performed in the factory method
    }

    /**
     * Factory method for generating DifferentiableSum DifferentiableFunctions. Note that this method is NOT guaranteed to
     * return a DifferentiableSum: if the values passed in can be represented more simply within this library, then they will be. The
     * only guarantee on the return type is that it will be some kind of DifferentiableFunction, and that it will be the
     * simplest available representation of the input.
     * @param terms The terms of the sum
     * @return A DifferentiableFunction representing the sum of the input terms
     */
    public static DifferentiableFunction construct(DifferentiableFunction... terms) {
        return construct(1, Arrays.asList(terms));
    }

    /**
     * Factory method for generating DifferentiableSum DifferentiableFunctions. Note that this method is NOT guaranteed to
     * return a DifferentiableSum: if the values passed in can be represented more simply within this library, then they will be. The
     * only guarantee on the return type is that it will be some kind of DifferentiableFunction, and that it will be the
     * simplest available representation of the input.
     * @param coefficient The coefficient of the entire sum
     * @param terms The terms of the sum
     * @return A DifferentiableFunction representing the sum of the input terms, multiplied by the coefficient
     */
    public static DifferentiableFunction construct(double coefficient, DifferentiableFunction... terms) {
        return construct(coefficient, Arrays.asList(terms));
    }

    /**
     * Factory method for generating DifferentiableSum DifferentiableFunctions. Note that this method is NOT guaranteed to
     * return a DifferentiableSum: if the values passed in can be represented more simply within this library, then they will be. The
     * only guarantee on the return type is that it will be some kind of DifferentiableFunction, and that it will be the
     * simplest available representation of the input.
     * @param coefficient The coefficient of the entire sum
     * @param inputTerms The terms of the sum
     * @return A DifferentiableFunction representing the sum of the input terms, multiplied by the coefficient
     */
    public static DifferentiableFunction construct(double coefficient, List<DifferentiableFunction> inputTerms) {
        if (coefficient == 0) {
            return ZERO;
        }

        List<DifferentiableFunction> terms = new ArrayList<>(inputTerms);

        Map<Class<? extends DifferentiableFunction>, List<DifferentiableFunction>> aggregatedTerms = new HashMap<>();


        for (int i = 0; i<terms.size();) {
            DifferentiableFunction term = terms.get(i++);

            if (!ZERO.equals(term)) {
                term = term.multiply(coefficient);
                if (term instanceof DifferentiableSum sumTerm) {
                    terms.addAll(sumTerm.getTerms());
                } else {
                    aggregateTerm(aggregatedTerms, term);
                }
            }
        }

        List<DifferentiableFunction> nonZeroTerms = new ArrayList<>();
        aggregatedTerms.values().forEach(nonZeroTerms::addAll);

        if (nonZeroTerms.isEmpty()) {
            return ZERO;
        }
        if (nonZeroTerms.size() == 1) {
            return nonZeroTerms.get(0);
        }

        return new DifferentiableSum(nonZeroTerms);
    }

    /**
     * Given a Map keyed by DifferentiableFunction Type to a corresponding List of DifferentiableFunctions of that type,
     * and a term to add; cleanly and statically add the term to the Map. This means finding its Type in the map,
     * iterating through the corresponding list until a term is found that it adds to cleanly, removing that term from
     * the list, and then calling this function again on the same (modified) map to insert the result of the clean
     * addition. If no clean additions can be performed within the corresponding list, the new term is simply added to
     * the list without recursion.
     * @param aggregatedTerms A Map containing Lists, each of a particular Type of DifferentiableFunction. The Map is
     *                        keyed by corresponding Types
     * @param term The term to be inserted into the map
     * @see DifferentiableSum#construct(double, List)
     */
    private static void aggregateTerm(Map<Class<? extends DifferentiableFunction>, List<DifferentiableFunction>> aggregatedTerms, DifferentiableFunction term) {
        List<DifferentiableFunction> aggregatedTermsOfThisType =
                aggregatedTerms.computeIfAbsent(term.getClass(), k-> new ArrayList<>());
        for (int i = 0; i < aggregatedTermsOfThisType.size(); i++) {
            DifferentiableFunction extantTerm = aggregatedTermsOfThisType.get(i);
            Optional<DifferentiableFunction> sum = extantTerm.addCleanly(term);
            if (sum.isPresent()) {
                aggregatedTermsOfThisType.remove(extantTerm);
                aggregateTerm(aggregatedTerms, sum.orElseThrow());
                return;
            }
        }
        aggregatedTermsOfThisType.add(term);
    }


    /**
     * @return an immutable copy of the terms of this sum
     */
    public List<DifferentiableFunction> getTerms() {
        return List.copyOf(terms);
    }

    @Override
    public DifferentiableFunction factorise() {
        List<DifferentiableFunction> commonFactors = new ArrayList<>(terms.get(0).getFactors());
        int i = 1;
        while (i < terms.size() && ! commonFactors.isEmpty()) {
            commonFactors.retainAll(terms.get(i).getFactors());
            i++;
        }
        List<DifferentiableFunction> firstRoundFactorisedTerms = new ArrayList<>();
        for (DifferentiableFunction unfactorisedTerm : terms) {
            DifferentiableFunction factorisedTerm = unfactorisedTerm;
            for (DifferentiableFunction commonFactor : commonFactors) {
                // This should never throw, since we're only dividing by factors that the terms themselves declared that they had
                factorisedTerm = factorisedTerm.divideCleanly(commonFactor).orElseThrow();
            }
            firstRoundFactorisedTerms.add(factorisedTerm);
        }
        // TODO further rounds of factorisation
        commonFactors.add(construct(1, firstRoundFactorisedTerms));
        return DifferentiableProduct.construct(1, commonFactors);
    }

    @Override
    protected DifferentiableFunction differentiateExplicitly() {
        List<DifferentiableFunction> differentiatedTerms = new ArrayList<>();

        for (DifferentiableFunction term : terms) {
            differentiatedTerms.add(term.differentiate());
        }

        return construct(getCoefficient(), differentiatedTerms);
    }

    @Override
    public Optional<DifferentiableFunction> divideCleanly(DifferentiableFunction divisor) {

        List<DifferentiableFunction> newTerms = new ArrayList<>();
        for (DifferentiableFunction term : terms) {
            Optional<DifferentiableFunction> dividedTerm = term.divide(divisor);
            if (dividedTerm.isEmpty()) {
                return Optional.empty();
            } else {
                newTerms.add(dividedTerm.orElseThrow());
            }
        }
        return Optional.of(construct(getCoefficient(), newTerms));
    }

    @Override
    public Optional<DifferentiableFunction> multiplyCleanly(DifferentiableFunction multiplicand) {
        List<DifferentiableFunction> multipliedTerms = new ArrayList<>();
        for (DifferentiableFunction term : terms) {
            Optional<DifferentiableFunction> multipliedTerm = term.multiplyCleanly(multiplicand);
            if (multipliedTerm.isEmpty()) {
                return Optional.empty();
            } else {
                multipliedTerms.add(multipliedTerm.orElseThrow());
            }
        }
        return Optional.of(construct(getCoefficient(), multipliedTerms));
    }

    @Override

    public @NotNull DifferentiableFunction multiply(DifferentiableFunction multiplicand) {
        List<DifferentiableFunction> multipliedTerms = new ArrayList<>();
        for (DifferentiableFunction term : terms) {
            multipliedTerms.add(term.multiply(multiplicand));
        }
        return construct(getCoefficient(), multipliedTerms);
    }

    @Override
    public Optional<DifferentiableFunction> subtractCleanly(DifferentiableFunction subtrahend) {
        for (DifferentiableFunction term : terms) {
            Optional<DifferentiableFunction> subtractedTerm = term.subtractCleanly(subtrahend);
            if (subtractedTerm.isPresent()) {
                List<DifferentiableFunction> updatedTerms = new ArrayList<>(terms);
                updatedTerms.remove(term);
                updatedTerms.add(subtractedTerm.orElseThrow());
                return Optional.of(construct(getCoefficient(), updatedTerms));
            }
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> addCleanly(DifferentiableFunction addend) {
        for (DifferentiableFunction term : terms) {
            Optional<DifferentiableFunction> addedTerm = term.addCleanly(addend);
            if (addedTerm.isPresent()) {
                List<DifferentiableFunction> updatedTerms = new ArrayList<>(terms);
                updatedTerms.remove(term);
                updatedTerms.add(addedTerm.orElseThrow());
                return Optional.of(construct(getCoefficient(), updatedTerms));
            }
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> exponentiateCleanly(double exponent) {
        return Optional.empty();
    }

    @Override
    public @NotNull String toString(int significantFigures) {
        StringBuilder string = new StringBuilder();

        string.append(terms.get(0).toString());
        for (DifferentiableFunction term : terms.subList(1, terms.size())) {
            string.append(" + %s".formatted(term.toString(significantFigures)));
        }
        return getCoefficient()==1? string.toString() : "%s(%s)".formatted(stripTrailingZeros(roundToSignificantFigures(getCoefficient(),significantFigures)), string.toString());
    }

    @Override
    public @NotNull DifferentiableFunction withCoefficient(double coefficient) {
        if (coefficient == getCoefficient()) {
            return this;
        }
        return construct(coefficient, terms);
    }


    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        DifferentiableSum differentiableSum = (DifferentiableSum) o;

        if (Double.compare(differentiableSum.getCoefficient(), getCoefficient()) != 0) {
            return false;
        }

        List<DifferentiableFunction> otherTerms = new ArrayList<>(differentiableSum.terms);
        if (otherTerms.size() != terms.size()) {
            return false;
        }
        for (DifferentiableFunction term : terms) {
            if (!otherTerms.remove(term)) {
                return false;
            }
        }
        return otherTerms.isEmpty();
    }

    @Override
    public int hashCode() {
        List<Integer> hashCodeList = new ArrayList<>();
        for (DifferentiableFunction term : terms) {
            hashCodeList.add(term.hashCode());
        }
        Collections.sort(hashCodeList);

        return Objects.hash("+", getCoefficient(), hashCodeList);
    }

    @Override
    public DifferentiableFunction withInputFunction(DifferentiableFunction newInputFunction) {
        return compose(newInputFunction);
    }

    @Override
    public DifferentiableFunction compose(DifferentiableFunction innerFunction) {
        List<DifferentiableFunction> composedTerms = new ArrayList<>(terms.size());
        for (DifferentiableFunction term : terms) {
            composedTerms.add(term.compose(innerFunction));
        }
        return construct(getCoefficient(), composedTerms);
    }

    @Override
    public List<DifferentiableFunction> commonFactors(List<DifferentiableFunction> potentialCommonFactors) {
        return null;
    }

    @Override
    public double applyAsDouble(double input) {
        double sum = 0;
        for (DifferentiableFunction term : terms) {
            sum += term.applyAsDouble(input);
        }
        return getCoefficient() * sum;
    }


}
