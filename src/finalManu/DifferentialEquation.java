package finalManu;

/**
 * Created by Manuel on 22/02/14.
 */
public abstract class DifferentialEquation {

    /*Class Differential equation. Represents a Differential equation with its derivatives.
      Every method returns the result of the evaluation of the differential equation or a
      derivative in a point x.
      Use and implementation of differential equations is needed for using the rungeKutta
      and Taylor methods in class Numeric.
      Generally, the form of a differential equation is y' = f(x,y). If that is the case, evaluate()
      would be referred to that equation, and the second derivative would be y'' = f'(x,y). The "first"
      derivative of the differential equation is actually the second derivative of the function y.
      @author: Manuel Acuña
      @author: Nicolás Rudolph
    */

    public DifferentialEquation(){};

    /*evaluate: evaluates the differential equation in a point x and in an initial condition y.
     *example = Math.cos(x) - Math.sin(y) + Math.pow(x,2);
     *@param x - point where the function is going to be evaluated.
      @param y - initial condition y0
     */

    abstract double evaluate(double x, double y);

    /*derivative: evaluates the differential equation´s derivative in a point -x- and in an initial condition -y-.
      @param x - point where the function is going to be evaluated.
      @param y - initial condition y0
      */

    abstract double derivative(double x, double y);

    /*secondDerivative: evaluates the differential equation's second derivative in a point -x- and in an initial condition -y-
    * @param x - point where the function is going to be evaluated.
      @param y - initial condition y0*/

    abstract double secondDerivative(double x, double y);

    /*thirdDerivatve: evaluates the differential equation's third derivative in a point -x- and in an initial condition -y-
    * @param x - point where the function is going to be evaluated.
      @param y - initial condition y0*/

    abstract double thirdDerivative(double x, double y);

    /*fourthDerivatve: evaluates the differential equation's fourth derivative in a point -x- and in an initial condition -y-
    * @param x - point where the function is going to be evaluated.
      @param y - initial condition y0*/

    abstract double fourthDerivative(double x, double y);


}
