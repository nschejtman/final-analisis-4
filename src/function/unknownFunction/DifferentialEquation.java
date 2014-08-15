package function.unknownFunction;

/**
 * Differential equation with two parameters x and y.
 */
public abstract class DifferentialEquation {
    /**
     * Evaluates the function in a (x, y) point.
     *
     * @param x x-coordinate
     * @param y y-coordinate
     * @return evaluation result
     */
    public abstract double evaluate(double x, double y);

    /**
     * Evaluates the function's derivative in a (x, y) point.
     *
     * @param x x-coordinate
     * @param y y-coordinate
     * @return evaluation result
     */
    public abstract double derivative(double x, double y);

    /**
     * Evaluates the function's 2nd order derivative in a (x, y) point.
     *
     * @param x x-coordinate
     * @param y y-coordinate
     * @return evaluation result
     */
    public abstract double secondDerivative(double x, double y);

    /**
     * Evaluates the function's 3rd order derivative in a (x, y) point.
     *
     * @param x x-coordinate
     * @param y y-coordinate
     * @return evaluation result
     */
    public abstract double thirdDerivative(double x, double y);

    /**
     * Evaluates the function's 4th order derivative in a (x, y) point.
     *
     * @param x x-coordinate
     * @param y y-coordinate
     * @return evaluation result
     */
    public abstract double fourthDerivative(double x, double y);


}
