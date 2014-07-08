package finalManu;

/**
 * Created by Manuel on 17/02/14.
 */
public abstract class Function {

    /*Class Function. Represents a Function with its derivatives. Every method returns
      the result of the evaluation of function or a derivative in a point x.
      @param x - point where the function is going to be evaluated.
      @author: Manuel Acuña
      @author: Nicolás Rudolph
    */

    public Function() {} ;

    /*evaluate: evaluates the function in a point. Implementation needed for all methods of the class Numeric
     * example = Math.pow(x,2) - 2;
     * @param x - point where the function is going to be evaluated.*/

    abstract double evaluate(double x);

    /* derivative: evaluates the derivative in a point. Implementation necessary
    * for: Newton-Raphson,bisection, integration(trapecios, simpson,romberg)
    * example = 2*x;
    * @param x - point where the function is going to be evaluated.
    *  */

    abstract double derivative(double x);

    /* secondDerivative: evaluates the second derivative in a point. Implementation necessary
    * for: Newton-Raphson and integration by trapecios;
    * example = 2;
    * @param x - point where the function is going to be evaluated.
    *  */

    abstract double secondDerivative(double x);

    /*thirdDerivative: evaluates the third derivative in a point. needed for solving differential equations
    * @param x - point where the function is going to be evaluated.*/

    abstract double thirdDerivative(double x);

    /*fourthDerivative: evaluates the fourth derivative in a point. needed for solving differential equations
      @param x - point where the function is going to be evaluated.
    * */

    abstract double fourthDerivative(double x);

    /* More derivatives can be implemented if necessary.  */
}
