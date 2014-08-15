package util;

import function.Function;

public abstract class Summation {


    /**
     * The sum of trapezoids returns the value resulting from trapezoidz'formula
     *
     * @param function the function to be integrated
     * @param a        lower limit of the integral
     * @param b        Upper limit of the integral
     * @param h        increment value
     */
    public static double trapeziusSum(Function function, double a, double b, double h) {
        double sum = 0;
        for (int i = 1; i < ((b - a) / h) - 1; i++)
            sum += function.evaluate(a + (i * h));
        return sum;
    }

    /**
     * The sum of trapezoids returns the value resulting from trapezoidz'formula
     *
     * @param f the function to be integrated
     * @param a lower limit of the integral
     * @param h increment value
     * @param n number of partitions, you must be an even number
     * @return
     */
    public static double trapeziusSum2(Function f, double a, double h, double n) {
        double sum = 0;
        for (int i = 1; i < n; i++)
            sum += f.evaluate(a + (i * h));
        return sum;
    }

    /**
     * Method that returns the sum of the even slots for simpson method
     *
     * @param function the function to be integrated
     * @param a        lower limit of the integral
     * @param b        Upper limit of the integral
     * @param n        number of partitions, you must be an even number
     * @return
     */
    public static double simpsonEven(Function function, double a, double b, double n) {
        double sum = 0;
        for (int i = 1; i <= (n / 2); i++)
            sum += function.evaluate(a + (2 * i * ((b - a) / n)));
        return sum;
    }

    /**
     * Method that returns the sum of the odd intervals for simpson method
     *
     * @param function the function to be integrated
     * @param a        lower limit of the integral
     * @param b        Upper limit of the integral
     * @param n        number of partitions, you must be an even number
     * @return
     */
    public static double simpsonOdd(Function function, double a, double b, double n) {
        double sum = 0;
        for (int i = 1; i <= (n / 2) - 1; i++)
            sum = sum + function.evaluate(a + (((2 * i) - 1) * ((b - a) / n)));
        return sum;
    }

}
