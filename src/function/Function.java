package function;

/**
 * Interface for an f(x) type function (only one variable). Has three methods: derivative & integral return a Function
 * class which is the result of derivative or integral of the given function. Evaluate returns a double corresponding to
 * the f(x) value for a given x value.
 */
public interface Function {

    /**
     * First order derivative for the function.
     * @return first derivative.
     */
    public Function derivative();

    /**
     * Integral for the function. To calculate the area simply use the formula:
     * Area = Function.integral(upper bound) - Function.integral(lower bound)
     * @return integral.
     */
    public Function integral();

    /**
     * Evaluates f(x) for a given x value.
     * @param x point to be evaluated.
     * @return result of the evaluation.
     */
    public double evaluate(double x);
}
