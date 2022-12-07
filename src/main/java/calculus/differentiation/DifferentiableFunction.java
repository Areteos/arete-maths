package calculus.differentiation;

import calculus.differentiation.functionTypes.Constant;
import calculus.differentiation.functionTypes.Monomial;
import net.jcip.annotations.Immutable;
import org.jetbrains.annotations.NotNull;

import java.util.*;
import java.util.function.DoubleUnaryOperator;

// TODO: NaturalLog, Sine, Cosine, and Tangent subclasses
/**
 * <p>An abstract class representing the family of functions of one variable that can be analytically differentiated.
 * The class provides a variety of methods for analytic function manipulation, construction, algebra, and
 * differentiation. This is noticeably faster and more accurate (especially for repeated differentiation) than numerical
 * methods.</p>
 * <p>Every DifferentiableFunction has a coefficient and an input function. This is because these concepts are universal
 * to all functions, and so can be implemented at a high level.</p>
 * <p>E.g. the function f(x)=x<sup>2</sup> is a Monomial with coefficient of 1 and no input function. By comparison,
 * g(x)=3(x+1)<sup>2</sup> is the same Monomial with the input function h(x)=x+1 and a coefficient of 3 </p>
 * <p>More specific mathematical components like exponents are less universal, and therefore handled at a concrete
 * implementation level.</p>
 * <p>All object instantiation of DifferentiableFunctions is intended to be done either through factory methods or
 *  subclass constructors where the constructed object is guaranteed to be the most mathematically simple interpretation
 *  of the input. E.g. there should be no way to construct the function t(x)=0x<sup>2</sup> as a
 *  DifferentiableFunction, attempts to do so should always result in the return of a Constant of value zero.</p>
 */
@Immutable
public abstract class DifferentiableFunction implements DoubleUnaryOperator {
    public static final int DEFAULT_SIGNIFICANT_FIGURES = 3;
    private final DifferentiableFunction inputFunction;
    private final double coefficient;
    private final static Map<DifferentiableFunction, DifferentiableFunction> differentiations = new HashMap<>();

    /**
     * The canonical constructor of all DifferentiableFunctions. Should never be called outside the constructors of
     * child classes.
     * @param inputFunction The inputFunction of this DifferentiableFunction
     * @param coefficient The coefficient of this DifferentiableFunction
     */
    protected DifferentiableFunction(DifferentiableFunction inputFunction, double coefficient) {
        if (Double.isInfinite(coefficient)) {
            throw new ArithmeticException("Infinite coefficient");
        }
        this.inputFunction = Objects.equals(inputFunction, Monomial.X)?null:inputFunction;
        this.coefficient = coefficient;
    }

    /**
     * <p>Find the first derivative of this function and return it. Note that the result may not necessarily be of the
     * same class as the original function. For example, the Monomial f(x)=x, when differentiated, returns the Constant
     * g(x)=1</p>
     * <p>Note that differentiation is evaluated lazily. Disregarding the coefficient, if an identical function has been
     * differentiated before, the previous result is returned. Each unique function will therefore only be explicitly
     * differentiated once per runtime.</p>
     * @return The derivative of this function
     */
    public final DifferentiableFunction differentiate() {
        DifferentiableFunction result = retrieveDifferentiationResult();
        if (result == null) {
            result = differentiateExplicitly();
            // If the dividing the result by this function's coefficient is undefined, that means our coefficient is
            // zero. In which case, differentiating this function should just give the zero function.
            storeDifferentiationResult(result.divide(coefficient).orElse(Constant.ZERO));
        } else {
            result = result.multiply(coefficient);
        }
        return result;
    }

    /**
     * Differentiate this function n times, and return the result. Note that this method cannot perform
     * antidifferentiation: values of n of 0 or less will simply return this function.
     * @param n The number of times to differentiate the function
     * @return The nth derivative of the function
     * @see DifferentiableFunction#differentiate() differentiate
     */
    public final DifferentiableFunction differentiate(int n) {
        DifferentiableFunction result = this;
        if (n < 0) {
            System.err.printf("WARNING: Negative input to repeat differentiation function: %s%n", n);
        }
        for (int i = 0; i < n; i++) {
            result = result.differentiate();
        }
        return result;
    }

    /**
     * Explicitly differentiates this function, without making use of lazy utilities.
     * @return The first derivative of this function
     * @see DifferentiableFunction#differentiate()
     */
    protected abstract DifferentiableFunction differentiateExplicitly();

    /**
     * Factorises this function and returns the result. This is currently exclusively useful for optimisation, and only
     * for some function types. Note that the result is exactly mathematically equivalent to the original function, but
     * may contain fewer computationally expensive calls to functions like exp().
     * @return A factorised version of this function. May simply return itself if no factorisation is possible.
     */
    public DifferentiableFunction factorise() {
        return this;
    }

    /**
     * Algebraically divide this function by another differentiable function. Will only return non-empty if the
     * operation is defined and the quotient is at least as simple as this function or the divisor.
     * @param divisor The function by which to divide this one
     * @return The quotient of the division. Will be empty if 0 is used as a divisor or the result is not clean.
     */
    public abstract Optional<DifferentiableFunction> divideCleanly(DifferentiableFunction divisor);


    /**
     * Algebraically divide this function by another differentiable function.
     * @param divisor The function by which to divide this one
     * @return The quotient of the division. Will be empty if 0 is used as a divisor.
     */
    public final Optional<DifferentiableFunction> divide(DifferentiableFunction divisor) {
        if (divisor instanceof Constant) {
            return divide(divisor.getCoefficient());
        }
        return Optional.of(
                divideCleanly(divisor).orElse(
                multiply(Monomial.construct(1, -1, divisor))));

    }

    /**
     * Algebraically divide this function by a number.
     * @param divisor The number by which to divide this function
     * @return The quotient of the division. Will be empty if 0 is used as a divisor.
     */
    public Optional<DifferentiableFunction> divide(double divisor) {
        if (divisor == 0) {
            return Optional.empty();
        }
        return Optional.of(withCoefficient(coefficient / divisor));
    }

    /**
     * Algebraically multiply this function by another differentiable function. Will return empty unless the product is
     * at least as simple as this function or the multiplicand.
     * @param multiplicand The function by which to multiply this one
     * @return The product of the multiplication, or empty if the result is not clean
     */
    public abstract Optional<DifferentiableFunction> multiplyCleanly(DifferentiableFunction multiplicand);

    /**
     * Algebraically multiply this function by another differentiable function.
     * @param multiplicand The function by which to multiply this one
     * @return The product of the multiplication
     */
    @NotNull
    public DifferentiableFunction multiply(DifferentiableFunction multiplicand) {
        if (multiplicand instanceof Constant) {
            return multiply(multiplicand.getCoefficient());
        }
        return multiplyCleanly(multiplicand).orElse(
                DifferentiableProduct.construct(this, multiplicand));
    }

    /**
     * Algebraically multiply this function by a number.
     * @param multiplicand The number by which to multiply this one
     * @return The product of the multiplication
     */
    @NotNull
    public final DifferentiableFunction multiply(double multiplicand) {
        return withCoefficient(coefficient * multiplicand);
    }

    /**
     * Algebraically subtract another differentiable function from this one. Will return empty unless the difference is
     * at least as simple as this function or the subtrahend.
     * @param subtrahend The function to subtract from this one
     * @return The product of the subtraction, or empty if the result is not clean
     */
    public abstract Optional<DifferentiableFunction> subtractCleanly(DifferentiableFunction subtrahend);

    /**
     * Algebraically subtract another differentiable function from this function.
     * @param subtrahend The function to subtract from this one
     * @return The difference of this subtraction
     */
    @NotNull
    public final DifferentiableFunction subtract(DifferentiableFunction subtrahend) {
        return subtractCleanly(subtrahend).orElse(
                DifferentiableSum.construct(this, subtrahend.multiply(-1)));
    }

    /**
     * Algebraically add another differentiable function to this one. Will return empty unless the sum is
     * at least as simple as this function or the addend.
     * @param addend The function to add to this one
     * @return The product of the addition, or empty if the result is not clean
     */
    public abstract Optional<DifferentiableFunction> addCleanly(DifferentiableFunction addend);

    /**
     * Algebraically add another differentiable function to this function.
     * @param addend The function to add to this one
     * @return The sum of this addition
     */
    @NotNull
    public final DifferentiableFunction add(DifferentiableFunction addend) {
        return addCleanly(addend).orElse(DifferentiableSum.construct(this, addend));
    }

    /**
     * Algebraically raise this function to an integer exponent.
     * @param exponent the exponent to which to raise this function
     * @return The result of the exponentiation. Will be empty if the exponent is negative and this function is zero
     */
    public final Optional<DifferentiableFunction> exponentiate(int exponent) {
        Optional<DifferentiableFunction> cleanDoubleResult = exponentiateCleanly((double) exponent);
        if (cleanDoubleResult.isPresent()) {
            return cleanDoubleResult;
        }
        DifferentiableFunction product = new Constant(1);
        if (exponent >= 0) {
            if (this.coefficient == 0) {
                return Optional.of(this);
            }
            for (int i = 0; i < exponent; i++) {
                product = product.multiply(this);
            }
            return Optional.of(product);
        } else {
            if (this.coefficient == 0) {
                return Optional.empty();
            }
            return Optional.of(Monomial.construct(1, -1, this.exponentiate(Math.abs(exponent)).orElseThrow()));
        }
    }

    /**
     * Algebraically raise this function to an integer exponent. Will return empty unless the power is at least as simple
     * as this function
     * @param exponent the exponent to which to raise this function
     * @return The result of the exponentiation. Will be empty if the exponent is negative and this function is zero
     */
    public final Optional<DifferentiableFunction> exponentiateCleanly(int exponent) {
        Optional<DifferentiableFunction> cleanDoubleResult = exponentiateCleanly((double) exponent);
        if (cleanDoubleResult.isPresent()) {
            return cleanDoubleResult;
        }
        DifferentiableFunction product = new Constant(1);
        if (exponent >= 0) {
            if (this.coefficient == 0) {
                return Optional.of(this);
            }
            for (int i = 0; i < exponent; i++) {
                Optional<DifferentiableFunction> result = product.multiplyCleanly(this);
                if (result.isEmpty()) {
                    return Optional.empty();
                } else {
                    product = result.orElseThrow();
                }
            }
            return Optional.of(product);
        } else {  // TODO investigate clean negative integer exponentiation
            return Optional.empty();
        }
    }

    /**
     * Algebraically raise this function to a rational exponent.
     * @param exponent the exponent to which to raise this function
     * @return The result of the exponentiation. Will be empty if the exponent is negative and this function is zero
     */
    @NotNull
    public final Optional<DifferentiableFunction> exponentiate(double exponent) {
        Optional<DifferentiableFunction> cleanDoubleResult = exponentiateCleanly(exponent);
        if (cleanDoubleResult.isPresent()) {
            return cleanDoubleResult;
        }
        if (exponent % 1 == 0) {
            return exponentiate((int) exponent);
        }
        if (exponent < 0 && coefficient==0) {
            return Optional.empty();
        }
        return Optional.of(Monomial.construct(1, exponent, this));
    }

    /**
     * Algebraically raise this function to a rational exponent. Will return empty unless the power is at least as simple
     * as this function
     * @param exponent the exponent to which to raise this function
     * @return The result of the exponentiation. Will be empty if the exponent is negative and this function is zero
     */
    public abstract Optional<DifferentiableFunction> exponentiateCleanly(double exponent);

    /**
     * Returns a String representation of this DifferentiableFunction, with all constant values rounded to the number
     * of significant figures given by DEFAULT_SIGNIFICANT_FIGURES
     * @return A String representation of the Differentiable Function
     * @see DifferentiableFunction#toString(int)
     * @see DifferentiableFunction#DEFAULT_SIGNIFICANT_FIGURES
     */
    @Override
    public final String toString() {
        return toString(DEFAULT_SIGNIFICANT_FIGURES);
    }

    /**
     * Returns a String representation of this DifferentiableFunction, with all constant values rounded to the given
     * number of significant figures (note that that is NOT the same thing as decimal places).
     * @param significantFigures Number of significant figures to round to
     * @return A string representation of this number, rounded to the given number of significant figures
     * @see DifferentiableFunction#toString()
     */
    @NotNull
    public abstract String toString(int significantFigures);

    /**
     * @return String representation of this function's input function
     * @param significantFigures Number of significant figures to round to
     * @see DifferentiableFunction#toString(int)
     */
    @NotNull
    public String inputFunctionToString(int significantFigures) {
        return inputFunction==null?"x":"(%s)".formatted(inputFunction.toString(significantFigures));
    }

    public double getCoefficient() {
        return coefficient;
    }

    /**
     * Returns a new DifferentiableFunction that is otherwise exactly the same as the present one, but with the given
     * coefficient.
     * @param coefficient The coefficient for the new function
     * @return A function identical to this one, but with the given coefficient
     */
    @NotNull
    public abstract DifferentiableFunction withCoefficient(double coefficient);

    @Override
    public abstract boolean equals(Object other);

    @Override
    public abstract int hashCode();

    /**
     * @return The exact value of this function's inputFunction field. Will be null if inputFunction is null
     */
    public DifferentiableFunction getInputFunction() {
        return inputFunction;
    }

    /**
     * @return An accurate representation of this function's inputFunction field that is guaranteed to be non-null.
     * Will return a Monomial representing the function f(x)=x if inputFunction is null.
     */
    public DifferentiableFunction getOperableInputFunction() {
        return Objects.requireNonNullElse(inputFunction, Monomial.X);
    }

    /**
     * Composes this function with another, returning a new composite function. The composite function will apply the
     * given function first on any input, then this one on the result to give its final result. The composite function
     * will also be a fully-fledged DifferentiableFunction, capable of being further used in the same manner as any
     * other. Properly handles recursion to ensure the nested inputFunction structures are respected.
     * @param innerFunction The inner function for the composition
     * @return A composite function of this function composed with the given inner function
     */
    public DifferentiableFunction compose(DifferentiableFunction innerFunction) {
        if (inputFunction == null) {
            return withInputFunction(innerFunction);
        } else {
            return withInputFunction(inputFunction.compose(innerFunction));
        }
    }

    /**
     * Creates a new differentiable function, identical to this one but with the given input function.
     * @param newInputFunction The input function to set
     * @return An identical function to this one, except with a different input function
     */
    public abstract DifferentiableFunction withInputFunction(DifferentiableFunction newInputFunction);

    /**
     * A thread-safe method to access the stored result of differentiating this function.
     * @return The prior result of differentiating this function
     */
    private synchronized DifferentiableFunction retrieveDifferentiationResult() {
        return differentiations.get(withCoefficient(1));
    }

    /**
     * A thread-safe method to store the result of differentiating this function.
     * @param result The result of differentiating this function
     */
    private synchronized void storeDifferentiationResult(DifferentiableFunction result) {
        differentiations.put(withCoefficient(1), result);
    }

    /**
     * Evaluates and returns a list of all factors of this function. Useful for factorisation and simplification. Every
     * factor in the list will be given with a coefficient of 1, and the coefficient of this function will NOT be
     * considered.
     * @return The factors of this function
     */
    public List<DifferentiableFunction> getFactors() {
        return List.of(this.withCoefficient(1));
    }

    /**
     * Given a list of factors, return a list of those which are factors in this function.
     * @param potentialCommonFactors Candidates for common factors between this DifferentiableFunction and another
     * @return A list of factors that are actually present in this DifferentiableFunction. This may include, for example, reducing the exponent of a monomial in the original list to find commonality
     */
    public abstract List<DifferentiableFunction> commonFactors(List<DifferentiableFunction> potentialCommonFactors);
}
