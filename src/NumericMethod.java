import function.Function;
import function.unknownFunction.DifferentialEquation;
import util.ExtendedMath;
import util.Summation;

/**
 * Class that implements numeric methods for root calculation (Newton-Raphson & Bisection), integral calculus
 * (Simpson, Trapeziums & Romberg) and differential equation evaluation (Taylor, Euler, 2nd order Runge-Kutta,
 * 4th order Runge-Kutta, Adams-Bashforth & Adams-Moulton)
 *
 * @author Diego Fabra
 * @author Nicolas Schejtman
 */
public abstract class NumericMethod {

    /**
     * Bisection numeric method to approximate the root of a f(x) function in a given interval. This method has three
     * conditions:
     * 1. There should be only one root in the interval
     * 2. f(from) * f(to) < 0 -> has a root in the interval
     * 3. df/dx != 0 for from<x<to
     * <p/>
     * It calculates the root by taking the mid value of the interval and redefining it as the new upper or lower bound
     * for a new interval, until the mid value is equal to the root.
     *
     * @param function function to be analyzed
     * @param from     lower bound
     * @param to       upper bound
     * @return root
     */
    public static double bisection(Function function, double from, double to, double epsilon) {
        //check df/dx != 0
        Function derivative = function.derivative();
        boolean isConstant = true;
        for (int i = 0; i < 500; i++) {
            if (derivative.evaluate(i) != 0) isConstant = false;
        }
        if (isConstant) throw new IllegalArgumentException("knownFunctions.Function is constant");

        //check f(from)*f(to) < 0
        final double yFrom = function.evaluate(from);
        final double yTo = function.evaluate(to);
        if (yFrom * yTo >= 0) throw new IllegalArgumentException("None or multiple roots in the interval");

        //Apply method
        final double middle = (from + to) / 2;
        final double yMiddle = function.evaluate(middle);

        if (Math.abs(yMiddle) < epsilon) return middle;
        else {
            if (yFrom * yMiddle < 0) return bisection(function, from, middle, epsilon);
            else return bisection(function, middle, to, epsilon);
        }
    }

    /**
     * Newton Raphson numeric method to approximate the root of a f(x) function in a given interval. This method has three
     * conditions:
     * 1. There should be only one root in the interval
     * 2. f(from) * f(to) < 0 -> has a root in the interval
     * 3. df/dx != 0 -> is not constant
     * <p/>
     * It calculates the root by calculating the x-intercept of either the tangent of the upper or lower bound of the
     * interval (which ever of the two that is within the interval), and redefining such x coordinate as the new upper
     * or lower bound of the interval, until the x-intercept is equal to the root.
     *
     * @param function function to be analyzed
     * @param from     lower bound
     * @param to       upper bound
     * @return root
     */
    public static double newtonRaphson(Function function, double from, double to, double epsilon) {
        //check df/dx != 0
        Function derivative = function.derivative();
        boolean isConstant = true;
        double i = from;
        while (i < to) {
            i = i + epsilon;
            if (derivative.evaluate(i) != 0) isConstant = false;
        }


        if (isConstant) throw new IllegalArgumentException("knownFunctions.Function is constant");

        //check f(from)*f(to) < 0
        final double yFrom = function.evaluate(from);
        final double yTo = function.evaluate(to);
        if (yFrom * yTo >= 0) throw new IllegalArgumentException("None or multiple roots in the interval");

        //Apply method
        double x0 = to;
        double x1 = x0 - (function.evaluate(x0) / function.derivative().evaluate(x0));
        if (x1 > to || x1 < from) {
            x0 = from;
            x1 = x0 - (function.evaluate(x0) / function.derivative().evaluate(x0));
        }

        return newtonRaphson(function, x1, epsilon);

    }

    private static double newtonRaphson(Function function, double x0, double epsilon) {
        double x1 = x0 - (function.evaluate(x0) / function.derivative().evaluate(x0));
        final double yMiddle = function.evaluate(x1);
        if (Math.abs(yMiddle) < epsilon) return x1;
        else return newtonRaphson(function, x1, epsilon);

    }

    /**
     * Taylor method for evaluating differential equations. Given a starting point (x0, yo), this method applies
     * Taylor's polynomial for function approximation to calculate the next point (x1, y1) in which x1 = x0 + h and
     * y1 = f(x0 + h) (where f is approximated with Taylor's polynomial), which will serve as the starting point for the
     * next iteration. The algorithm iterates until the next x value is the value given for evaluation and returns the
     * corresponding y value.
     *
     * @param diff differential equation
     * @param x0   starting x coordinate
     * @param y0   starting y coordinate
     * @param h    increase for x value in each iteration
     * @param x    value to be evaluated
     * @return y value for x
     */
    public static double taylorDiff(DifferentialEquation diff, double x0, double y0, double h, double x) {
        while (x0 < x) {
            //Calculate the next point
            double y1 =
                    y0 +
                            h * diff.evaluate(x0, y0) +
                            ((java.lang.Math.pow(h, 2)) / ExtendedMath.factorial(2)) * diff.derivative(x0, y0) +
                            ((java.lang.Math.pow(h, 3)) / ExtendedMath.factorial(3)) * diff.secondDerivative(x0, y0) +
                            ((java.lang.Math.pow(h, 4)) / ExtendedMath.factorial(4)) * diff.thirdDerivative(x0, y0) +
                            ((java.lang.Math.pow(h, 5)) / ExtendedMath.factorial(5)) * diff.fourthDerivative(x0, y0);
            double x1 = x0 + h;

            //Redefine the starting point
            y0 = y1;
            x0 = x1;
        }
        return y0;
    }

    /**
     * Euler method for evaluating differential equations. Given a starting point (x0, yo), this method applies Euler's
     * iteration to calculate the next point (x1, y1) which will then serve as the starting point for the next iteration.
     * The iteration stops when the starting x-coordinate is equal to the x value to be evaluated and returns the
     * corresponding y-coordinate.
     *
     * @param diff differential equation
     * @param x0   starting x coordinate
     * @param y0   starting y coordinate
     * @param h    increase for x value in each iteration
     * @param x    value to be evaluated
     * @return y value for x
     */
    public static double euler(DifferentialEquation diff, double x0, double y0, double h, double x) {
        while (x0 < x) {
            //Calculate the next point
            double y1 = y0 + h * diff.evaluate(x0, y0); //Euler iteration
            double x1 = x0 + h;

            //Redefine the starting point
            x0 = x1;
            y0 = y1;
        }
        return y0;
    }

    /**
     * Runge-Kutta method for evaluating differential equations. Given a starting point (x0, yo), this method applies
     * 2nd order Runge-Kutta iteration to calculate the next point (x1, y1) which will then serve as the starting point
     * for the next iteration. The iteration stops when the starting x-coordinate is equal to the x value to be
     * evaluated and returns the corresponding y-coordinate.
     *
     * @param diff differential equation
     * @param x0   starting x coordinate
     * @param y0   starting y coordinate
     * @param h    increase for x value in each iteration
     * @param x    value to be evaluated
     * @return y value for x
     */
    public static double rungeKuttaO2(DifferentialEquation diff, double x0, double y0, double h, double x) {
        while (x0 < x) {
            //Calculate the next point
            double f1 = h * diff.evaluate(x0, y0);
            double f2 = h * diff.evaluate(x0 + h, y0 + f1);

            double y1 = y0 + 0.5 * (f1 + f2); //2nd Order Runge-Kutta iteration
            double x1 = x0 + h;

            //Redefine the starting point
            x0 = x1;
            y0 = y1;
        }
        return y0;
    }

    /**
     * Runge-Kutta method for evaluating differential equations. Given a starting point (x0, yo), this method applies
     * 4th order Runge-Kutta iteration to calculate the next point (x1, y1) which will then serve as the starting point
     * for the next iteration. The iteration stops when the starting x-coordinate is equal to the x value to be
     * evaluated and returns the corresponding y-coordinate.
     *
     * @param diff differential equation
     * @param x0   starting x coordinate
     * @param y0   starting y coordinate
     * @param h    increase for x value in each iteration
     * @param x    value to be evaluated
     * @return y value for x
     */
    public static double rungeKuttaO4(DifferentialEquation diff, double x0, double y0, double h, double x) {
        while (x0 < x) {
            //Calculate the next point
            double f1 = h * diff.evaluate(x0, y0);
            double f2 = h * diff.evaluate(x0 + 0.5 * h, y0 + 0.5 * f1);
            double f3 = h * diff.evaluate(x0 + 0.5 * h, y0 + 0.5 * f2);
            double f4 = h * diff.evaluate(x0 + h, y0 + f3);

            double y1 = y0 + (1/6.) * (f1 + 2 * f2 + 2 * f3 + f4); //4th Order Runge-Kutta iteration
            double x1 = x0 + h;

            //Redefine the starting point
            x0 = x1;
            y0 = y1;
        }
        return y0;
    }

    /**
     * Adams-Bashforth method for evaluating differential equations. Given a starting point (x1, y1) and its predecesor
     * (x0, yo), this method applies Adams-Bashforth iteration to calculate the next point (x2, y2) which will then serve
     * as the starting point for the next iteration. The iteration stops when the starting x-coordinate is equal to the
     * x value to be evaluated and returns the corresponding y-coordinate.
     *
     * @param diff differential equation
     * @param x0   predecessor x-coordinate
     * @param y0   predecessor y-coordinate
     * @param x1   starting x-coordinate
     * @param y1   starting y-coordinate
     * @param h    increase for x value in each iteration
     * @param x    value to be evaluated
     * @return y value for x
     */
    public static double adamsBashforth(DifferentialEquation diff, double x0, double y0, double x1, double y1, double h, double x) {
        while (x1 < x) {
            //Calculate the next point
            double f0 = diff.evaluate(x0, y0);
            double f1 = diff.evaluate(x1, y1);

            double y2 = y1 + h * (1.5 * f1 - 0.5 * f0); //Adams-Bashforth iteration
            double x2 = x1 + h;

            //Redefine the starting point and its predecessor
            y0 = y1;
            x0 = x1;
            y1 = y2;
            x1 = x2;
        }
        return y1;
    }

    /**
     * Adams-Moulton method for evaluating differential equations. Given a starting point (x1, y1) and its predecesor
     * (x0, yo), this method applies Adams-Bashforth iteration to calculate the next point (x2, y2) and then
     * Adams-Moulton iteration to reduce Adams-Bashforth's calculation error. This point will then serve as the starting
     * point for the next iteration. The iteration stops when the starting x-coordinate is equal to the x value to be
     * evaluated and returns the corresponding y-coordinate.
     *
     * @param diff differential equation
     * @param x0   predecessor x-coordinate
     * @param y0   predecessor y-coordinate
     * @param x1   starting x-coordinate
     * @param y1   starting y-coordinate
     * @param h    increase for x value in each iteration
     * @param x    value to be evaluated
     * @return y value for x
     */
    public static double adamsMoulton(DifferentialEquation diff, double x0, double y0, double x1, double y1, double h, double x) {
        while (x1 < x) {
            //Calculate the next point
            double x2 = x1 + h;
            double y2 = adamsBashforth(diff, x0, y0, x1, y1, h, x);

            double f1 = diff.evaluate(x1, y1);
            double f2 = diff.evaluate(x2, y2);

            y2 = y1 + h * (0.5 * f1 + 0.5 * f2);

            //Redefine the starting point and its predecessor
            y0 = y1;
            x0 = x1;
            y1 = y2;
            x1 = x2;
        }
        return y1;
    }


    /**
     * The composite trapezoidal rule or trapezoidal rule is a way of approximating a definite integral using n trapezoids.
     * In the formulation of this method assumes that f is continuous and positive on the interval [a, b].
     *
     * @param function the function to be integrated
     * @param from     lower limit of the integral
     * @param to       Upper limit of the integral
     * @param h        increment value
     * @return
     */
    public static double trapeziums(Function function, double from, double to, double h) {
        double n = (to - from) / h;
        double root = (h / 2) * (function.evaluate(from) + function.evaluate(to) + (2 * Summation.trapeziusSum(function, from, to, h)));
        // double error = ((to - from) / 12) * Math.pow(h, 2) * maxOfSecondDerivative(function, from, to);
        //System.out.print("answer" + root );
        return root;
    }

    /**
     * In numerical integration, a way of approximating a definite integral for an interval [a, b] is using the trapezoidal rule, that is,
     * that on each subinterval in which divides [a, b] approximates f by a polynomial first grade, and then calculate the integral as the sum of the areas of the trapezoids formed in these subintervals.
     * The method used for Simpson's rule follows the same philosophy, but approaching the subintervals of f by quadratic polynomials.
     *
     * @param function the function to be integrated
     * @param from     lower limit of the integral
     * @param to       Upper limit of the integral
     * @param h        increment value
     * @return value of the integral
     */

    public static double simpson(Function function, double from, double to, double h) {


        double n = (to - from) / h;
        double result = (h / 3) * (function.evaluate(from) + 2 * Summation.simpsonOdd(function, from, to, n) + 4 * Summation.simpsonEven(function, from, to, n) + function.evaluate(to));
        //double error = ((to - from) / 180) * Math.pow(h, 4) * maxOfFourthDerivative(function, from, to);
        //System.out.println("result: " + result);
        //System.out.println("error: " + " +/- " + error);
        return result;
    }

    /**
     * @param function the function to be integrated
     * @param from     lower limit of the integral
     * @param to       Upper limit of the integral
     * @param row      double that set the row of the last romberg term, the term that is going to be returned
     * @param column   double that set the column of the last romberg term, the term that is going to be returned
     * @return last romberg term
     */
    public static double integrateByRomberg(Function function, double from, double to, int row, int column) { // row and column are "coordinates" that set what romberg term is going to be returned.
        int k = 0;     //represents the number corresponding to the h (column);
        double hk = (to - from) / (Math.pow(2, k));
        double n = (to - from) / hk;
        double[][] matrix = new double[row + 1][];  // this matrix is going to be filled with romberg terms
        int rowCapacity = 1;  //i use this variable to build the matrix
        for (int j = matrix.length - 1; j >= 0; j--) { // this for loop builds the matrix
            matrix[j] = new double[column + rowCapacity];  //the last position is going to have size =  column + rowCapacity, the row bfore the last one: column + rowCapacity+1, and so on
            rowCapacity++;
        }
        for (int l = 0; l < matrix[0].length; l++) { // this for loop fills the first row of the matrix with the trapezium terms that are needed for getting romberg terms
            double t0k = (hk / 2) * (function.evaluate(from) + function.evaluate(to) + 2 * Summation.trapeziusSum2(function, from, hk, n));   //Pasar h y no N
            matrix[0][l] = t0k;
            k++;
            hk = (to - from) / (Math.pow(2, k));
            n = (to - from) / hk;
        }
        int position; // this integer moves along the row in which the method inserts the romberg terms
        double romberg = 0;
        for (int d = 0; d < matrix.length - 1; d++) { // d represents the previous row to the one in which the method inserts the romberg terms
            position = 0;   // each time the program enters the loop, position returns to zero, in order to go across the row again.
            for (int g = 1; g < matrix[d].length; g++) {
                romberg = (Math.pow(4, d + 1) * matrix[d][g] - matrix[d][g - 1]) / ((Math.pow(4, d + 1)) - 1);
                matrix[d + 1][position] = romberg;
                position++;
            }
        }
//        double error = ((Math.pow(to - from, 3)) / (12 * Math.pow(row, 2))) * maxOfSecondDerivative(function, from, to);
        //System.out.println("error: " + " +/- " + error);
        return matrix[row][column];  // returns the romberg term given by the coordinates passed through parameters
    }

    /**
     * Method that calculates the max of the second derivative in a given interval.
     * Needed for calculating the error in trapezium and Romberg integration methods.
     *
     * @param function - Function to be used, for finding approximations to its roots.
     * @param from     - double that represents the first number to be evaluated, the initial point.
     * @param to       - double that represents the end of the interval.
     *                 return max.
     */

    private static double maxOfSecondDerivative(Function function, double from, double to) {
        double max = function.derivative().evaluate(from);
        for (double i = from; i <= to; i += 0.001) {
            double result = function.derivative().evaluate(i);
            if (result > max) {
                max = result;
            }
        }
        return max;
    }

    /**
     * Method that calculates the max of the fourth derivative in a given interval.
     * Needed for calculating the error in Simpson's integration method.
     *
     * @param function - Function to be used, for finding approximations to its roots.
     * @param from     - double that represents the first number to be evaluated, the initial point.
     * @param to       - double that represents the end of the interval.
     *                 return max.
     */
    private static double maxOfFourthDerivative(Function function, double from, double to) {

        final Function der4 = function.derivative().derivative().derivative().derivative();
        if (der4 != null) {

            double max = der4.evaluate(from);
            for (double i = from; i <= to; i += 0.001) {
                double result = der4.evaluate(i);
                if (result > max) {
                    max = result;
                }
            }
            return max;
        } else {
            return 0;
        }
    }
}
