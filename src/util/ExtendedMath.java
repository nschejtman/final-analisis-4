package util;

/**
 * Extended mathematical functions (factorial)
 */
public abstract class ExtendedMath {

    /**
     * Returns the factorial of a number
     *
     * @param number
     * @return factorial
     */
    public static double factorial(double number) {
        if (number == 0) return 1;
        else return number * factorial(number - 1);
    }

}
