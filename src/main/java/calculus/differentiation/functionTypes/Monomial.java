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

import calculus.differentiation.DifferentiableFunction;
import calculus.differentiation.DifferentiableSum;
import calculus.differentiation.DifferentiableProduct;
import functions.MathsFunctions;
import types.tuples.Pair;
import net.jcip.annotations.Immutable;
import org.jetbrains.annotations.NotNull;


import java.util.*;

import static functions.MathsFunctions.roundToSignificantFigures;
import static functions.MiscFunctions.stripTrailingZeros;


/**
 * <p>An immutable class representing a differentiable function of the form:</p>
 * &emsp&emsp f(x) = ax<sup>b</sup>
 * <p>where a and b are constant.</p>
 * Expansions of sums are considered simpler than exponentiated sums: ax<sup>2</sup>+axy+ay<sup>2</sup> is considered
 * simpler than a(x+y)<sup>2</sup>
 */
@Immutable
public class Monomial extends DifferentiableFunction {
    final double exponent;

    /**
     * A static representation of the function f(x) = x, very useful in a variety of places.
     */
    public final static Monomial X = new Monomial(1, 1, null);

    /**
     * Constructs a Monomial with exactly the information given, no further processing performed. Can only be called
     * safely from within the factory method of the Monomial class.
     * @param coefficient The coefficient for this Monomial
     * @param exponent The exponent of this Monomial
     * @param inputFunction The inputFunction of this Monomial, or null if x
     * @see Monomial#construct(double, double, DifferentiableFunction)
     */
    private Monomial(double coefficient, double exponent, DifferentiableFunction inputFunction) {
        super(inputFunction, coefficient);
        this.exponent = exponent;
    }


    /**
     * <p>
     * Factory method for generating Monomial DifferentiableFunctions. Note that this method is NOT guaranteed to return
     * a Monomial: if the values passed in can be represented more simply within this library, then they will be. The
     * only guarantee on the return type is that it will be some kind of DifferentiableFunction, and that it will be the
     * simplest available representation of the input. Expansions of sums are considered simpler than exponentiated
     * sums: ax<sup>2</sup>+axy+ay<sup>2</sup> is considered simpler than a(x+y)<sup>2</sup></p>
     * <p>Since the Monomial function is simply exponentiating some input by a constant and then multiplying the power
     * by a constant, there are many potential inputs to this method that would result in a non-Monomial output</p>
     * @param coefficient The coefficient of the intended Monomial
     * @param exponent The exponent of the intended Monomial
     * @param inputFunction The input function to the intended Monomial, or null if x
     * @return A DifferentiableFunction representing the function f(x)=a(g(x))<sup>b</sup>, where a is the coefficient,
     * b is the exponent, and g(x) is the input function
     */
    @NotNull
    public static DifferentiableFunction construct(double coefficient, double exponent, DifferentiableFunction inputFunction) {
        if (coefficient == 0) {  // 0 times anything is 0
            return Constant.ZERO;
        }
        if (exponent == 0) {  // Any input function raised to 0 is 1
            return new Constant(coefficient);
        }
        if (inputFunction != null) {
            if (exponent == 1) {  // If the exponent is 1, this is really just a multiple of the input function
                return inputFunction.multiply(coefficient);
            }
            if (inputFunction instanceof Constant) {  // If you always input a constant, you always output a constant
                if (inputFunction.getCoefficient() < 0 && exponent % 1 != 0) {
                    throw new ArithmeticException("Complex number operations not supported");
                } if (inputFunction.getCoefficient() == 0 && exponent < 0) {
                    throw new ArithmeticException("Division by zero");
                }
                return new Constant(coefficient * Math.pow(inputFunction.getCoefficient(), exponent));
            } else if (inputFunction instanceof Monomial inputMonomial) {
                return construct(  // A monomial with a monomial input can be multiplied out into a simpler monomial
                        coefficient * Math.pow(inputMonomial.getCoefficient(), exponent),
                        exponent * inputMonomial.exponent,
                        inputMonomial.getInputFunction()
                );
            } else if (inputFunction instanceof NaturalExponent) {
                return NaturalExponent.construct(  // Natural exponents can always be further exponentiated neatly
                        coefficient * Math.pow(inputFunction.getCoefficient(), exponent),
                        inputFunction.getOperableInputFunction().withCoefficient(inputFunction.getOperableInputFunction().getCoefficient() * exponent)
                );
            } else if (inputFunction instanceof DifferentiableProduct) {  // Simplicity is debatable here, but bringing
                // the exponent into an exponentiated product should result in a lower maximum nesting level of functions overall
                List<DifferentiableFunction> currentFactors = inputFunction.getFactors();
                List<DifferentiableFunction> newFactors = new ArrayList<>(currentFactors.size());
                for (DifferentiableFunction currentFactor : currentFactors) {
                    newFactors.add(construct(1, exponent, currentFactor));
                }
                return DifferentiableProduct.construct(inputFunction.getCoefficient() * coefficient, newFactors);
            } else if (exponent % 1 == 0 && exponent > 0 && inputFunction instanceof DifferentiableSum) {
                // If we have a positive integer exponent and a sum as the input, simply expand the sum.
                // Since we check for positivity, and have already checked for a constant (i.e. zero) input, this
                // exponentiation should never return empty.
                return inputFunction.exponentiate((int) exponent).orElseThrow().multiply(coefficient);
            }

            // If all of the above passes, bring the input function's coefficient out to the main function's
            if (inputFunction.getCoefficient() != 1) {
                coefficient *= Math.pow(inputFunction.getCoefficient(), exponent);
                inputFunction = inputFunction.withCoefficient(1);
            }
        }

        // And, finally, construct a new monomial and return it
        return new Monomial(coefficient, exponent, inputFunction);
    }


    @Override
    protected DifferentiableFunction differentiateExplicitly() {
        double newCoefficient = getCoefficient() * exponent;
        double newExponent = exponent - 1;

        DifferentiableFunction newMainFunction = construct(newCoefficient, newExponent, getInputFunction());

        if (getInputFunction() == null) {
            return newMainFunction;
        } else {
            return DifferentiableProduct.construct(getInputFunction().differentiate(), newMainFunction);
        }
    }

    @Override
    public Optional<DifferentiableFunction> divideCleanly(DifferentiableFunction divisor) {

        if (divisor instanceof Constant) {  // Constant division is easy
            return super.divide(divisor.getCoefficient());
        } else if (divisor instanceof Monomial monomialDivisor) {  // Monomial division might be clean
            if (Objects.equals(getInputFunction(), monomialDivisor.getInputFunction())) {
                return Optional.of(construct(  // Very easy if input functions are equal
                        getCoefficient() / monomialDivisor.getCoefficient(),
                        exponent - monomialDivisor.exponent,
                        getInputFunction()));
            } else {  // Otherwise there is one more path we can take: Simply exponentiating both input functions and then dividing them together
                double outerExponent;
                int dividendInnerExponent;
                int divisorInnerExponent;
                if (exponent == monomialDivisor.exponent) {
                    dividendInnerExponent = 1;
                    divisorInnerExponent = 1;
                    outerExponent = exponent;
                } else if (exponent%1==0 && monomialDivisor.exponent%1==0) {
                    // If both exponents are integers we can easily simplify the problem
                    outerExponent = MathsFunctions.findGreatestCommonFactor((int) exponent, (int) monomialDivisor.exponent);
                    dividendInnerExponent = (int) (exponent / outerExponent);
                    divisorInnerExponent = (int) (monomialDivisor.exponent / outerExponent);
                } else {
                    // If not, we'll have to do some work to turn the exponents into integers, and our outer exponent
                    // will not be an integer.
                    // This should always be defined, since a monomial with exponent 0 should not exist
                    Pair<Integer, Integer> exponentFraction = MathsFunctions.simplifyFraction(exponent, monomialDivisor.exponent).orElseThrow();
                    dividendInnerExponent = exponentFraction.first();
                    divisorInnerExponent = exponentFraction.second();
                    outerExponent = exponent / dividendInnerExponent;
                }

                Optional<DifferentiableFunction> newInputFunction;
                try {
                    newInputFunction =
                            getOperableInputFunction().exponentiate(dividendInnerExponent).orElseThrow().divideCleanly(
                                    monomialDivisor.getOperableInputFunction().exponentiate(divisorInnerExponent).orElseThrow());
                } catch (NoSuchElementException e) {
                    return Optional.empty();
                }
                if (newInputFunction.isEmpty()) {
                    return Optional.empty();
                }
                return Optional.of(construct(
                        getCoefficient() / monomialDivisor.getCoefficient(),
                        outerExponent,
                        newInputFunction.orElseThrow()));
            }
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> multiplyCleanly(DifferentiableFunction multiplicand) {

        if (multiplicand instanceof Constant) {  // Multiplication with a constant is easy
            return Optional.of(super.multiply(multiplicand.getCoefficient()));
        } else if (multiplicand instanceof Monomial monomialMultiplicand) {  // Otherwise can only multiply cleanly with another Monomial
            if (Objects.equals(getInputFunction(), monomialMultiplicand.getInputFunction())) {
                return Optional.of(construct(  // If input functions are equal, then multiplication is easy
                        getCoefficient() * multiplicand.getCoefficient(),
                        exponent + ((Monomial) multiplicand).exponent,
                        getInputFunction()
                ));

            } else {  // Otherwise we'll have to work with the input functions directly
                double outerExponent;
                int multiplierInnerExponent, multiplicandInnerExponent;
                if (exponent == monomialMultiplicand.exponent) {
                    multiplierInnerExponent = 1;
                    multiplicandInnerExponent = 1;
                    outerExponent = exponent;
                } else if (exponent % 1 == 0 && monomialMultiplicand.exponent % 1 == 0) {
                    outerExponent = MathsFunctions.findGreatestCommonFactor((int) exponent, (int) monomialMultiplicand.exponent);
                    multiplierInnerExponent = (int) (exponent / outerExponent);
                    multiplicandInnerExponent = (int) (monomialMultiplicand.exponent / outerExponent);
                } else {
                    Pair<Integer, Integer> exponentFraction = MathsFunctions.simplifyFraction(exponent, monomialMultiplicand.exponent).orElseThrow();
                    multiplierInnerExponent = exponentFraction.first();
                    multiplicandInnerExponent = exponentFraction.second();
                    outerExponent = exponent / multiplierInnerExponent;
                }

                Optional<DifferentiableFunction> newInputFunction;
                try {
                    newInputFunction =
                            getOperableInputFunction().exponentiate(multiplierInnerExponent).orElseThrow().multiplyCleanly(
                                    monomialMultiplicand.getOperableInputFunction().exponentiate(multiplicandInnerExponent).orElseThrow());
                } catch (NoSuchElementException e) {
                    return Optional.empty();
                }
                if (newInputFunction.isEmpty()) {
                    return Optional.empty();
                }
                return Optional.of(construct(
                        getCoefficient() * monomialMultiplicand.getCoefficient(),
                        outerExponent,
                        newInputFunction.orElseThrow()));


            }
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> subtractCleanly(DifferentiableFunction subtrahend) {

        if (
                subtrahend instanceof Monomial monomialSubtrahend &&
                exponent == monomialSubtrahend.exponent &&
                Objects.equals(getInputFunction(), subtrahend.getInputFunction())) {

            return Optional.of(construct(getCoefficient() - subtrahend.getCoefficient(), exponent, getInputFunction()));
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> addCleanly(DifferentiableFunction addend) {
        if (
                addend instanceof Monomial monomialAddend &&
                exponent == monomialAddend.exponent &&
                Objects.equals(getInputFunction(), addend.getInputFunction())) {

            return Optional.of(construct(getCoefficient() + addend.getCoefficient(), exponent, getInputFunction()));
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> exponentiateCleanly(double exponent) {
        return Optional.of(construct(Math.pow(getCoefficient(), exponent), this.exponent * exponent, getInputFunction()));
    }

    @Override
    public @NotNull String toString(int significantFigures) {
        return "%s%s".formatted(
                getCoefficient()!=1? stripTrailingZeros(roundToSignificantFigures(getCoefficient(),significantFigures)):"",
                inputFunctionToString(significantFigures))
                + (exponent!=1?"^%s".formatted(stripTrailingZeros(roundToSignificantFigures(exponent, significantFigures))):"");
    }

    @Override
    public @NotNull DifferentiableFunction withCoefficient(double coefficient) {
        if (coefficient == this.getCoefficient()) {
            return this;
        }
        return construct(coefficient, exponent, getInputFunction());
    }

    @Override
    public double applyAsDouble(double input) {
        if (getInputFunction() != null) {
            input = getInputFunction().applyAsDouble(input);
        }
        if (input == 0 && exponent < 0) {
            throw new ArithmeticException("Base of zero with negative exponent");
        }
        if (input < 0 && exponent % 1 != 0) {
            throw new ArithmeticException("Complex number operations not supported");
        }

        return getCoefficient() * Math.pow(input, exponent);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Monomial monomial = (Monomial) o;
        return Double.compare(monomial.getCoefficient(), getCoefficient()) == 0 && Double.compare(monomial.exponent, exponent) == 0 && Objects.equals(getInputFunction(), monomial.getInputFunction());
    }

    @Override
    public int hashCode() {
        return Objects.hash(getCoefficient(), exponent, Objects.hashCode(getInputFunction()));
    }

    @Override
    public DifferentiableFunction withInputFunction(DifferentiableFunction newInputFunction) {
        return construct(getCoefficient(), exponent, newInputFunction);
    }

    @Override
    public List<DifferentiableFunction> commonFactors(List<DifferentiableFunction> potentialCommonFactors) {
        Double greatestCommonExponent = null;
        for (DifferentiableFunction potentialCommonFactor : potentialCommonFactors) {
            if (potentialCommonFactor.getClass() == Monomial.class) {
                if (Objects.equals(getInputFunction(), potentialCommonFactor.getInputFunction())) {
                    if (greatestCommonExponent == null) {
                        greatestCommonExponent = exponent;
                    }
                    double potentialExponent = ((Monomial) potentialCommonFactor).exponent;
                    if (potentialExponent < greatestCommonExponent) {
                        greatestCommonExponent = potentialExponent;
                    }
                }
            }
        }
        if (greatestCommonExponent == null) {
            return List.of();
        } else {
            return List.of(construct(1, greatestCommonExponent, getInputFunction()));
        }
    }

}
