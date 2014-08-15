package function.knownFunction;

import function.Function;

/**
 * Generic polynomial function in the form f(x) = sum(i=0, i=n) ai x^i, where ai are the coefficients and n is the order
 * of the function.
 */
public class PolynomialFunction implements Function {
    int order;
    double[] coefficients;

    /**
     * Constructor
     *
     * @param order
     */
    public PolynomialFunction(int order) {
        this.order = order;
        coefficients = new double[order + 1];
    }

    /**
     * Returns the order of a polynomial function. This is the maximum exponent of the function. For example, consider:
     * f(x) = x^3 + 4x^2 - 3x + 5
     * the order of this function is 3.
     *
     * @return order
     */
    public int getOrder() {
        return order;
    }

    /**
     * Returns the coefficient (ai) for a given exponent value. For example, consider:
     * f(x) = x^3 + 4x^2 - 3x + 5
     * the coefficient for 3 is 1; for 2, 4; for 1, -3; and for 0, 5
     *
     * @param exponent
     * @return coefficient
     */
    public double getCoefficient(int exponent) {
        return coefficients[exponent];
    }

    /**
     * Sets the coefficient (ai) for a given exponent value. For example, consider:
     * f(x) = x^3 + 4x^2 - 3x + 5
     * the coefficient for 3 is 1; for 2, 4; for 1, -3; and for 0, 5.
     *
     * @param exponent
     * @return coefficient
     */
    public void setCoefficient(int exponent, double value) {
        coefficients[exponent] = value;
    }

    /**
     * Returns the polynomial function's derivative: another polynomial function of order = n-1 (one less than the
     * order of the actual function).
     *
     * @return derivative
     */
    @Override
    public Function derivative() {
        PolynomialFunction deriv = new PolynomialFunction(order - 1);
        for (int i = 0; i <= order; i++) {
            deriv.setCoefficient(i, coefficients[i + 1]);
        }
        return deriv;
    }

    /**
     * Returns the polynomial function's derivative: another polynomial function of order = n+1 (one more than the
     * order of the actual function). For the a0 coefficient it generates a random value, so do not consider this value.
     *
     * @return derivative
     */
    @Override
    public Function integral() {
        PolynomialFunction integ = new PolynomialFunction(order + 1);
        integ.setCoefficient(0, Math.random());
        for (int i = 1; i <= order + 1; i++) {
            integ.setCoefficient(i, coefficients[i - 1]);
        }
        return integ;
    }

    /**
     * Evaluates the polynomial function in a x point.
     *
     * @param x point to be evaluated.
     * @return evaluation result
     */
    @Override
    public double evaluate(double x) {
        double y = 0;
        for (int i = 0; i <= order; i++) {
            y = y + coefficients[i] * Math.pow(x, i);
        }
        return y;
    }
}
