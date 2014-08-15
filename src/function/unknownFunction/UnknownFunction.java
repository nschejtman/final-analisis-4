package function.unknownFunction;

/**
 * Function determined by evaluation points from itself and its derivatives rather than a f(x) equation.
 */
public abstract class UnknownFunction {
    /**
     * Evaluates the function in a x-coordinate.
     *
     * @param x x-coordinate for evaluation
     * @return evaluation result
     */
    public abstract double evaluate(double x);

    /**
     * Evaluates the function's derivative in a x-coordinate.
     *
     * @param x x-coordinate for evaluation
     * @return evaluation result
     */
    public abstract double firstDerivative(double x);

    /**
     * Evaluates the function's 2nd order derivative in a x-coordinate.
     *
     * @param x x-coordinate for evaluation
     * @return evaluation result
     */
    public abstract double secondDerivative(double x);

    /**
     * Evaluates the function's 3rd order derivative in a x-coordinate.
     *
     * @param x x-coordinate for evaluation
     * @return evaluation result
     */
    public abstract double thirdDerivative(double x);

    /**
     * Evaluates the function's 4th order derivative in a x-coordinate.
     *
     * @param x x-coordinate for evaluation
     * @return evaluation result
     */
    public abstract double fourthDerivative(double x);
}
