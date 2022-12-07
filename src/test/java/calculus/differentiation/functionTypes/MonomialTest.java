package calculus.differentiation.functionTypes;

import calculus.differentiation.DifferentiableProduct;
import calculus.differentiation.DifferentiableSum;
import org.junit.jupiter.api.Test;
import calculus.differentiation.DifferentiableFunction;
import org.junit.jupiter.api.function.Executable;

import static org.junit.jupiter.api.Assertions.*;

class MonomialTest {
    DifferentiableFunction xSquared = Monomial.construct(1, 2, null);
    DifferentiableFunction xCubed = Monomial.construct(1, 3, null);
    DifferentiableFunction sqrtX = Monomial.construct(1, 0.5, null);
    DifferentiableFunction oneOverX = Monomial.construct(1, -1, null);
    DifferentiableFunction oneOverTwoSqrtX = Monomial.construct(0.5, -0.5, null);

    DifferentiableFunction sixX = Monomial.construct(6, 1, null);

    DifferentiableFunction threeXSquared = Monomial.construct(3, 2, null);
    DifferentiableFunction fourXCubed = Monomial.construct(4, 3, null);
    DifferentiableFunction fiveSqrtX = Monomial.construct(5, 0.5, null);
    DifferentiableFunction sixOverX = Monomial.construct(6, -1, null);
    DifferentiableFunction sevenOverSqrtX = Monomial.construct(7, -0.5, null);

    @Test
    void construct() {
        DifferentiableFunction xPlusOneSquared = Monomial.construct(1, 2, DifferentiableSum.construct(Monomial.X, new Constant(1)));
        DifferentiableFunction xSquaredPlusTwoXPlusOne = DifferentiableSum.construct(xSquared, Monomial.construct(2, 1, null), new Constant(1));
        assertEquals(xSquaredPlusTwoXPlusOne, xPlusOneSquared);

        DifferentiableFunction eToTheXAllSquared = Monomial.construct(1, 2, NaturalExponent.construct(1, null));
        DifferentiableFunction eToThe2X = NaturalExponent.construct(1, Monomial.construct(2, 1, null));
        assertEquals(eToThe2X, eToTheXAllSquared);

        DifferentiableFunction threeSquared = Monomial.construct(1, 2, new Constant(3));
        DifferentiableFunction nine = new Constant(9);
        assertEquals(nine, threeSquared);

        DifferentiableFunction fiveSqrtXAllCubed = Monomial.construct(1, 3, fiveSqrtX);
        DifferentiableFunction oneHundredAndTwnetyFiveXToThePowerOfThreeOverTwo = Monomial.construct(125, 1.5, null);
        assertEquals(fiveSqrtXAllCubed, oneHundredAndTwnetyFiveXToThePowerOfThreeOverTwo);

        DifferentiableFunction xTimesEToTheXAllSquared = Monomial.construct(1, 2, DifferentiableProduct.construct(Monomial.X, NaturalExponent.construct(1, null)));
        DifferentiableFunction xSquaredTimeEToTheTwoX = DifferentiableProduct.construct(Monomial.construct(1, 2, null), NaturalExponent.construct(1, Monomial.construct(2, 1, null)));
        assertEquals(xSquaredTimeEToTheTwoX, xTimesEToTheXAllSquared);
    }

    @Test
    void differentiate() {
        DifferentiableFunction derivativeOfXCubed = xCubed.differentiate();
        DifferentiableFunction secondDerivativeOfXCubed = derivativeOfXCubed.differentiate();

        DifferentiableFunction derivativeOfSqrtX = sqrtX.differentiate();

        for (double i = -10000; i < 10000; i+=0.5) {
            assertEquals(threeXSquared.applyAsDouble(i), derivativeOfXCubed.applyAsDouble(i));
            assertEquals(sixX.applyAsDouble(i), secondDerivativeOfXCubed.applyAsDouble(i));

            if (i <= 0) {
                double finalI = i;
                assertThrows(ArithmeticException.class, () -> oneOverTwoSqrtX.applyAsDouble(finalI));
            } else {
                assertEquals(oneOverTwoSqrtX.applyAsDouble(i), derivativeOfSqrtX.applyAsDouble(i));
            }
        }

        assertInstanceOf(Constant.class, secondDerivativeOfXCubed.differentiate());
    }

    @Test
    void apply() {
        for (double i = -10000; i < 10000; i+=0.5) {
            final double finalI = i;
            assertEquals(Math.pow(i, 2), xSquared.applyAsDouble(i), 1e-10);
            assertEquals(Math.pow(i, 3), xCubed.applyAsDouble(i), 1e-10);
            assertEquals(3 * Math.pow(i, 2), threeXSquared.applyAsDouble(i), 1e-10);
            assertEquals(4 * Math.pow(i, 3), fourXCubed.applyAsDouble(i), 1e-10);

            if (i < 0) {
                assertThrows(ArithmeticException.class, ()->sqrtX.applyAsDouble(finalI));
                assertThrows(ArithmeticException.class, ()->fiveSqrtX.applyAsDouble(finalI));
            } else {
                assertEquals(Math.sqrt(i), sqrtX.applyAsDouble(i), 1e-10);
                assertEquals(5 * Math.sqrt(i), fiveSqrtX.applyAsDouble(i), 1e-10);

            }

            if (i == 0) {
                assertThrows(ArithmeticException.class, ()->oneOverX.applyAsDouble(finalI));
                assertThrows(ArithmeticException.class, ()->sixOverX.applyAsDouble(finalI));
            } else {
                assertEquals(1/i, oneOverX.applyAsDouble(i), 1e-10);
                assertEquals(6/i, sixOverX.applyAsDouble(i), 1e-10);
            }

            if (i <= 0) {
                assertThrows(ArithmeticException.class, ()->oneOverTwoSqrtX.applyAsDouble(finalI));
                assertThrows(ArithmeticException.class, ()->sevenOverSqrtX.applyAsDouble(finalI));
            } else {
                assertEquals(0.5/Math.sqrt(i), oneOverTwoSqrtX.applyAsDouble(i), 1e-10);
                assertEquals(7/Math.sqrt(i), sevenOverSqrtX.applyAsDouble(i), 1e-10);
            }

        }
    }
}