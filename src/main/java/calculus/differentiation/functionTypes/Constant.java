package calculus.differentiation.functionTypes;

import calculus.differentiation.DifferentiableFunction;
import calculus.differentiation.DifferentiableSum;
import functions.MathsFunctions;
import functions.MiscFunctions;
import net.jcip.annotations.Immutable;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.Objects;
import java.util.Optional;

/**
 * <p>An immutable class representing a differentiable function of the form:<p>
 * &emsp&emsp f(x) = c
 * <p>where c is a constant.
 */
@Immutable
public final class Constant extends DifferentiableFunction {
    /**
     * A static representation of the Constant function f(x) = 0, very useful in a variety of places.
     */
    public static final Constant ZERO = new Constant(0);

    /**
     * Creates a new Constant DifferentiableFunction
     * @param value the value of this constant function
     * @see Constant
     */
    public Constant(double value) {
        super(null, value);
    }

    @Override
    protected DifferentiableFunction differentiateExplicitly() {
        return ZERO;
    }

//    @Override
//    public DifferentiableFunction divide(DifferentiableFunction divisor) {
//        Class<? extends DifferentiableFunction> divisorClass = divisor.getClass();
//        if (divisorClass == Constant.class) {
//            return new Constant(getCoefficient() / divisor.getCoefficient());
//        }
//        return null;
//    }


    @Override
    public Optional<DifferentiableFunction> divideCleanly(DifferentiableFunction divisor) {
        if (divisor instanceof Constant) {
            if (ZERO.equals(divisor)) {
                return Optional.empty();
            }
            return Optional.of(new Constant(getCoefficient() / divisor.getCoefficient()));
        } else if (divisor instanceof DifferentiableSum) {
            return Optional.empty();  // a/Σb is not as simple as a or Σb
        }
        return super.divide(divisor);  // Any other extant DifferentiableFunction divided into a constant should result in a simple change of coefficient and exponent
    }

    @Override
    public Optional<DifferentiableFunction> multiplyCleanly(DifferentiableFunction multiplicand) {
        return Optional.of(multiplicand.multiply(getCoefficient()));
    }

    @Override
    public Optional<DifferentiableFunction> subtractCleanly(DifferentiableFunction subtrahend) {
        if (subtrahend instanceof Constant) {  // TODO trig identities
            return Optional.of(new Constant(getCoefficient() - subtrahend.getCoefficient()));
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> addCleanly(DifferentiableFunction addend) {
        if (addend instanceof Constant) {  // TODO trig identities
            return Optional.of(new Constant(getCoefficient() + addend.getCoefficient()));
        }
        return Optional.empty();
    }

    @Override
    public Optional<DifferentiableFunction> exponentiateCleanly(double exponent) {
        return Optional.of(new Constant(Math.pow(getCoefficient(), exponent)));
    }

    @Override
    public @NotNull String toString(int significantFigures) {
        return MiscFunctions.stripTrailingZeros(MathsFunctions.roundToSignificantFigures(getCoefficient(), significantFigures));
    }

    @Override
    public @NotNull DifferentiableFunction withCoefficient(double coefficient) {
        return new Constant(coefficient);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Constant constant = (Constant) o;
        return Double.compare(constant.getCoefficient(), getCoefficient()) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(getCoefficient());
    }

    @Override
    public DifferentiableFunction compose(DifferentiableFunction innerFunction) {
        return this;
    }

    @Override
    public DifferentiableFunction withInputFunction(DifferentiableFunction newInputFunction) {
        return this;
    }

    @Override
    public List<DifferentiableFunction> commonFactors(List<DifferentiableFunction> potentialCommonFactors) {
        return List.of();
    }

    @Override
    public double applyAsDouble(double input) {
        return getCoefficient();
    }
}
