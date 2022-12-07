package calculus.differentiation.functionTypes;

import calculus.differentiation.DifferentiableFunction;
import calculus.differentiation.DifferentiableSum;
import calculus.differentiation.DifferentiableProduct;
import net.jcip.annotations.Immutable;
import org.jetbrains.annotations.NotNull;


import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.concurrent.CopyOnWriteArrayList;

import static calculus.differentiation.functionTypes.Constant.ZERO;
import static functions.MathsFunctions.roundToSignificantFigures;
import static functions.MiscFunctions.stripTrailingZeros;


/**
 * <p>An immutable class representing a differentiable function of the form:</p>
 * &emsp&emsp f(x) = ae<sup>x</sup>
 * <p>where a is constant and e is Euler's number.</p>
 * <p>ae<sup>x+y</sup> will always be considered simpler than be<sup>x</sup>e<sup>y</sup></p>
 */
@Immutable
public class NaturalExponent extends DifferentiableFunction {

    /**
     * Constructs a NaturalExponent with exactly the information given, no further processing performed. Can only be called
     * safely from within the factory method of the NaturalExponent class.
     * @param coefficient The coefficient for this NaturalExponent
     * @param inputFunction The inputFunction of this NaturalExponent, or null if x
     * @see NaturalExponent#construct(double, DifferentiableFunction)
     */
    private NaturalExponent(double coefficient, DifferentiableFunction inputFunction) {
        super(inputFunction, coefficient);
    }

    /**
     * <p>
     * Factory method for generating NaturalExponent DifferentiableFunctions. Note that this method is NOT guaranteed to
     * return a NaturalExponent: if the values passed in can be represented more simply within this library, then they will be. The
     * only guarantee on the return type is that it will be some kind of DifferentiableFunction, and that it will be the
     * simplest available representation of the input. <p>ae<sup>x + y</sup> will always be considered simpler than
     * be<sup>x</sup>e<sup>y</sup></p></p>
     * @param coefficient The coefficient of the intended NaturalExponent
     * @param inputFunction The input function to the intended NaturalExponent, or null if x
     * @return A DifferentiableFunction representing the function f(x)=ae<sup>g(x)</sup>, where a is the coefficient,
     * e is Euler's number, and g(x) is the input function
     */
    @NotNull
    public static DifferentiableFunction construct(double coefficient, DifferentiableFunction inputFunction) {
        if (coefficient == 0) {
            return ZERO;
        }
        if (inputFunction instanceof Constant) {
            return new Constant(coefficient * Math.exp(inputFunction.getCoefficient()));
        }

        if (inputFunction instanceof DifferentiableSum inputSum) {   // Given a sum as the input function...
            List<DifferentiableFunction> retainedTerms = new ArrayList<>();
            for (DifferentiableFunction term : inputSum.getTerms()) {
                if (term instanceof Constant) {  // ...pare out all the Constant terms...
                    coefficient *= Math.exp(term.getCoefficient());  // ...and incorporate them into the coefficient
                } else {
                    retainedTerms.add(term);
                }
            }
            inputFunction = DifferentiableSum.construct(1, retainedTerms);
        }
        return new NaturalExponent(coefficient, inputFunction);
    }

    @Override
    protected DifferentiableFunction differentiateExplicitly() {
        if (getInputFunction() == null) {
            return this;
        } else {
            return DifferentiableProduct.construct(getInputFunction().differentiate(), this);
        }
    }

    @Override
    public Optional<DifferentiableFunction> divideCleanly(DifferentiableFunction divisor) {
        if (divisor instanceof Constant) {
            return super.divide(divisor.getCoefficient());
        } else if (divisor instanceof NaturalExponent) {
            double newCoefficient = getCoefficient() / divisor.getCoefficient();
            DifferentiableFunction newInputFunction = getOperableInputFunction().subtract(divisor.getOperableInputFunction());

            return Optional.of(construct(newCoefficient, newInputFunction));
        }
        return  Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> multiplyCleanly(DifferentiableFunction multiplicand) {

        if (multiplicand instanceof Constant) {
            return Optional.of(super.multiply(multiplicand.getCoefficient()));
        } else if (multiplicand instanceof NaturalExponent) {
            double newCoefficient = getCoefficient() * multiplicand.getCoefficient();
            DifferentiableFunction newInputFunction = getOperableInputFunction().add(multiplicand.getOperableInputFunction());

            return Optional.of(construct(newCoefficient, newInputFunction));
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> subtractCleanly(DifferentiableFunction subtrahend) {
        if (subtrahend instanceof NaturalExponent && Objects.equals(getInputFunction(), subtrahend.getInputFunction())) {
            return Optional.of(construct(getCoefficient() - subtrahend.getCoefficient(), getInputFunction()));
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> addCleanly(DifferentiableFunction addend) {
        if (addend instanceof NaturalExponent && Objects.equals(getInputFunction(), addend.getInputFunction())) {
            return Optional.of(construct(getCoefficient() + addend.getCoefficient(), getInputFunction()));
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> exponentiateCleanly(double exponent) {
        return Optional.of(construct(Math.pow(getCoefficient(), exponent), getOperableInputFunction().multiply(exponent)));
    }


    @Override
    public @NotNull String toString(int significantFigures) {
        return "%se^%s".formatted(
                getCoefficient()!=1? stripTrailingZeros(roundToSignificantFigures(getCoefficient(),significantFigures)):"",
                inputFunctionToString(significantFigures));
    }

    @Override
    public @NotNull DifferentiableFunction withCoefficient(double coefficient) {
        return construct(coefficient, getInputFunction());
    }

    @Override
    public double applyAsDouble(double input) {
        if (getInputFunction() != null) {
            input = getInputFunction().applyAsDouble(input);
        }

        return getCoefficient() * Math.exp(input);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        NaturalExponent naturalExponent = (NaturalExponent) o;
        return Double.compare(naturalExponent.getCoefficient(), getCoefficient()) == 0 && Objects.equals(getInputFunction(), naturalExponent.getInputFunction());
    }

    @Override
    public int hashCode() {
        return Objects.hash(getCoefficient(), Objects.hashCode(getInputFunction()));
    }

    @Override
    public DifferentiableFunction withInputFunction(DifferentiableFunction newInputFunction) {
        return construct(getCoefficient(), newInputFunction);
    }

    @Override
    public List<DifferentiableFunction> commonFactors(List<DifferentiableFunction> potentialCommonFactors) {

        DifferentiableFunction myInputFunction = getOperableInputFunction();
        DifferentiableFunction myInputFunctionWithoutCoefficient = myInputFunction.withCoefficient(1);

        Double greatestCommonExponentCoefficient = null;

        for (DifferentiableFunction potentialCommonFactor : potentialCommonFactors) {
            if (potentialCommonFactor.getClass() == NaturalExponent.class) {
                DifferentiableFunction otherInputFunction = potentialCommonFactor.getOperableInputFunction();

                if (Objects.equals(myInputFunctionWithoutCoefficient, otherInputFunction.withCoefficient(1))) {
                    if (greatestCommonExponentCoefficient == null) {
                        greatestCommonExponentCoefficient = myInputFunction.getCoefficient();
                    }
                    if (otherInputFunction.getCoefficient() < greatestCommonExponentCoefficient) {
                        greatestCommonExponentCoefficient = otherInputFunction.getCoefficient();
                    }
                }
            }
        }
        if (greatestCommonExponentCoefficient == null) {
            return List.of();
        } else {
            return List.of(construct(1, myInputFunction.withCoefficient(greatestCommonExponentCoefficient)));
        }
    }



    /**
     * Creates a DifferentiableFunction representing a Gaussian distribution with the given standard deviation and mean
     * @param standardDeviation Standard deviation of the distribution
     * @param mean Mean of the distribution
     * @return A Gaussian distribution as a DifferentiableFunction
     */
    public static DifferentiableFunction getGaussianDistribution(double standardDeviation, double mean) {
        double exponentCoefficient = - 0.5 / Math.pow(standardDeviation, 2);

        DifferentiableFunction innerFunction = DifferentiableSum.construct(exponentCoefficient,
                Monomial.construct(1, 2, null),
                Monomial.construct(-2 * mean, 1, null),
                new Constant(Math.pow(mean, 2)
        ));

        double coefficient = 1 / (standardDeviation * Math.sqrt(2 * Math.PI));

        return NaturalExponent.construct(coefficient, innerFunction);
    }
}
