package finalManu;

import java.util.ArrayList;
import java.util.Stack;

/**
 * Created by Manuel on 17/02/14.
 */
public class Numeric {


    private Numeric(){ }   // never instantiate this class

    /** Newton - Raphson. In numerical analysis, Newton–Raphson method, named
     * after Isaac Newton and Joseph Raphson, is a method for finding successively
     * better approximations to the roots (or zeroes) of a real-valued function.
     * the private method checkNRConditions, evaluates if the function passed through parameters
     * satisfies the newton-raphson conditions, which are:
     * f(from)*f(to) <0
     * f'(x) != 0 for every x belonging to the interval [from,to]
     * f has no inflexion points
     * if the function does not satisfy the conditions, the method breaks and returns an error;
     * executeNewtonRaphson is another private method, a recursive method that contains all the logic
     * of newton-raphson.
     * @param f - Function to be used, for finding aproximations to its roots.
     * @param from - double that represents the first number to be evaluated, the initial point.
     * @param to - double that represents the end of the interval.
     * @return root, or -1 if the function didn´t satisfy the conditions
     */

    public static double newtonRaphson(Function f, double from,double to){
       double root;
        if(checkNRConditions(f,from,to)){
           double x0 = -1;
           if(f.evaluate(from)*f.secondDerivative(from) > 0) x0 = from;
           if(f.evaluate(to)*f.secondDerivative(to) > 0) x0 = to;
           root = executeNewtonRaphson(f, x0);
           System.out.println("la raiz de la funcion ingresada, por el metodo de newton-raphson es: " + root);
           return root;
        }else{
            System.out.println(" La funcion no cumple las condiciones que requiere el metodo de Newton-Raphson."+"\n"+
            "No se pueden obtener las raices");
            return -1;
        }
    }


    /** Integration-Trapezium. In numerical analysis, the trapezium rule is a technique for approximating
     * a definite integral.
     * the private mehtod maxOfSecondDerivative, it is used for calculating the error of this method, the posible variation
     * that the result achieved by trapezium formula can have.
     * @param f - Function to be used, for finding aproximations to its roots.
     * @param a - double that represents the first number to be evaluated, the initial point.
     * @param b - double that represents the end of the interval.
     * @return the result of the integral./
     * */

    public static double integrateByTrapezium(Function f, double a, double b, double h){
        double n = (b-a)/h;
        double result = (h/2)*(f.evaluate(a)+f.evaluate(b)+2*Summatory.trapecios(f,a,h,n));   //trapezium formula
        double error = ((b-a)/12)*Math.pow(h,2)* maxOfSecondDerivative(f, a, b);
        System.out.println("result: " + result + " +/- " + error);
        return result;
    }

    /** Integration-Simpson. In numerical analysis, Simpson's rule is a method for numerical integration,
     * the numerical approximation of definite integrals.
     * @param f - Function to be integrated
     * @param a - double that represents the first integration limit.
     * @param b - double that represents the second integration limit.
     * @return result of the integral.*/

    public static double integrateBySimpson(Function f, double a, double b, double h){
        double n = (b-a)/h;
        double result = (h/3)*((f.evaluate(a)+2*Summatory.simpsonOdd(f,a,h,n)+4*Summatory.simpsonEven(f, a, h, n) + f.evaluate(b)));  //Simpson formula
        double error = ((b-a)/180)*Math.pow(h,4)*maxOfFourthDerivative(f,a,b);
        System.out.println("result: "+ result);
        System.out.println("error: " + " +/- "+ error);
        return result;
    }

    /** Integration-Romberg. Romberg's method  is used to estimate the definite integral
     * by applying Richardson extrapolation repeatedly on the trapezium rule or the
     * rectangle rule (midpoint rule). The estimates generate a triangular array.
     * The method is implemented by using a two-dimension array. In that way, the user can determine
     * the romberg term to be returned.
     * @param f - Function to be integrated
     * @param a - double that represents the first integration limit.
     * @param b - double that represents the second integration limit.
     * @param row - double that set the row of the last romberg term, the term that is going to be returned
     * @param column - double that set the column of the last romberg term, the term that is going to be returned
     * @return last romberg term
     * */

    public static double integrateByRomberg(Function f, double a, double b,int row, int column){ // row and column are "coordinates" that set what romberg term is going to be returned.
        int k = 0;     //represents the number corresponding to the h (column);
        double hk = (b-a)/(Math.pow(2,k));
        double n = (b-a)/hk;
        double[][] matrix = new double[row+1][];  // this matrix is going to be filled with romberg terms
        int rowCapacity = 1;  //i use this variable to build the matrix
        for(int j = matrix.length-1; j>=0 ;j--){ // this for loop builds the matrix
            matrix[j] = new double[column+rowCapacity];  //the last position is going to have size =  column + rowCapacity, the row bfore the last one: column + rowCapacity+1, and so on
            rowCapacity++;
        }
        for(int l = 0;l<matrix[0].length;l++){ // this for loop fills the first row of the matrix with the trapezium terms that are needed for getting romberg terms
            double t0k = (hk/2)*(f.evaluate(a) + f.evaluate(b) + 2*Summatory.trapecios(f,a,hk,n));
            matrix[0][l] = t0k;
            k++;
            hk = (b-a)/(Math.pow(2,k));
            n = (b-a)/hk;
        }
        int position; // this integer moves along the row in which the method inserts the romberg terms
        double romberg = 0;
        for(int d = 0;d<matrix.length-1;d++){ // d represents the previous row to the one in which the method inserts the romberg terms
            position = 0;   // each time the program enters the loop, position returns to zero, in order to go across the row again.
                for(int g = 1;g<matrix[d].length;g++){
                        romberg = (Math.pow(4,d+1)*matrix[d][g] - matrix[d][g-1])/((Math.pow(4,d+1))-1);
                        matrix[d+1][position] = romberg;
                        position++;
                }
        }
        double error = ((Math.pow(b-a,3))/(12*Math.pow(row,2)))*maxOfSecondDerivative(f,a,b);
        System.out.println("error: " + " +/- " + error);
        return matrix[row][column];  // returns the romberg term given by the coordinates passed through parameters
    }

    /** Differential Equation - Taylor. This method solves differential equations by aplying the Taylor
     * polinomium. The method returns a set of points that belong to the function y. The implementation of the four derivatives
     * of the Differential equation is necessary for a correct result.
     * @param yp - Differential Equation to be solved.
     * @param x0 - first value of the interval to be used
     * @param xn - last value of the interval to be used
     * @param h - growing interval of x0
     * @param y0 - initial value or condition, a point of the function y that is known beforehand
     * @return - ArrayList with the set of points*/


    public static ArrayList<Double> solveDiffEqByTaylor(DifferentialEquation yp, double x0, double xn, double h, double y0){
        ArrayList<Double> points = new ArrayList<Double>();
        while(x0<xn){
            double ykplus1 = y0 + h*yp.evaluate(x0, y0) + ((Math.pow(h,2))/factorial(2))*yp.derivative(x0, y0) +
                        ((Math.pow(h,3))/factorial(3))*yp.secondDerivative(x0, y0) + ((Math.pow(h,4))/factorial(4))*yp.thirdDerivative(x0,y0) +
                        ((Math.pow(h,5))/factorial(5))*yp.fourthDerivative(x0,y0);
            points.add(ykplus1);
            y0 = ykplus1;
            x0 = x0+h;
        }
        return points;
    }

    /** Differential Equation - Euler. First-order numerical procedure for solving ordinary differential
     * equations (ODEs) with a given initial value. The method returns a set of points that belong to the function y.
     * The formula for this method is: yk+1 = yk + h*(f(xk,yk). This method uses private method getF1(h,xk,yk); which returns
     * h*f(xk,yk);
     * @param yp - Differential Equation to be solved.
     * @param x0 - first value of the interval to be used
     * @param xn - last value of the interval to be used
     * @param h - growing interval of x0
     * @param y0 - initial value or condition, a point of the function y that is known beforehand
     * @return - ArrayList with the set of points*/

    public static ArrayList<Double> euler(DifferentialEquation yp, double x0, double xn, double h, double y0){
        ArrayList<Double> points = new ArrayList<Double>();
        while(x0<xn){
            double f1 = getF1(yp,h,x0,y0);
            double ykplus1 = y0 + h*f1;
            points.add(ykplus1);
            y0 = ykplus1;
            x0 = x0+h;
        }
        return points;
    }

    /** Differential Equations - Runge-Kutta2ndOrder. In numerical analysis, the Runge–Kutta methods are an important
     * family of implicit and explicit iterative methods, which are used in temporal discretization for the
     * approximation of solutions of ordinary differential equations. This is Runge-Kutta/2nd Order. Its formula is
     * yk+1 = yk +  h*(f(xk,yk) +  h*f(x+(1/2)*h,y+(1/2)*f1)), being f1 = h*(f(xk,yk). RungeKutta2ndOrder uses private method getF1 and
     * getF2. getF1(see euler method), getF2(h,xk,yk,f1) = h*f(x+(1/2)*h,y+(1/2)*f1)).
     * @param yp - Differential Equation to be solved.
     * @param x0 - first value of the interval to be used
     * @param xn - last value of the interval to be used
     * @param h - growing interval of x0
     * @param y0 - initial value or condition, a point of the function y that is known beforehand
     * @return - ArrayList with the set of points  */

    public static ArrayList<Double> rungeKutta2ndOrder(DifferentialEquation yp, double x0, double xn, double h, double y0){
        ArrayList<Double> points = new ArrayList<Double>();
        while(x0<xn){
            double f1 = getF1(yp, h, x0, y0);
            double f2 = getF2(yp, h, f1, x0, y0);
            double ykplus1 = y0 + ((f1+f2)/2);
            points.add(ykplus1);
            y0 = ykplus1;
            x0 = x0+h;
        }
        return points;
    }

    /** Differential Equations - Runge-Kutta4thOrder. In numerical analysis, the Runge–Kutta methods are an important
     * family of implicit and explicit iterative methods, which are used in temporal discretization for the
     * approximation of solutions of ordinary differential equations. This is Runge-Kutta/4th Order.
     * Its formula is: yk+1 = yk + (f1+2*f2+2*f3+f4). f1 given by the private method getF1(xk,yk) = h*(f(xk,yk) ; f2 given by
     * the private method giveF2(xk,yk,f1) =  h*f(x+(1/2)*h,y+(1/2)*f1)); f3 given by private method getF3(xk,yk,f2) = h*f(x+(1/2)*h,y+(1/2)*f2));
     * f4 given by the private method getF4(xk,yk,f3) = h*f(x+h,y+f3));
     * @param yp - Differential Equation to be solved.
     * @param x0 - first value of the interval to be used
     * @param xn - last value of the interval to be used
     * @param h - growing interval of x0
     * @param y0 - initial value or condition, a point of the function y that is known beforehand
     * @return - ArrayList with the set of points
     * */

    public static ArrayList<Double> rungeKutta4thOrder(DifferentialEquation yp, double x0, double xn, double h,double y0){
        ArrayList<Double> points = new ArrayList<Double>();
        while(x0<xn){
            double f1 = getF1(yp, h, x0, y0);
            double f2 = getF2(yp, h, f1, x0, y0);
            double f3 = getF3(yp, h, f2, x0, y0);
            double f4 = getF4(yp, h, f3, x0, y0);
            double ykplus1 = ((f1+2*f2+2*f3+f4)/6) + y0;
            points.add(ykplus1);
            y0 = ykplus1;
            x0 = x0+h;
        }
        return points;
    }

    public static int factorial(int a){
        if(a==1) return 1;
        else return a*factorial(a-1);

    }

    /** getF1. this method returns the f1 value needed for Runge-Kutta and euler methods
     * @param de - differential equation
     * @param h - growing interval of x
     * @param x - x that the DE uses
     * @param y - y initial condition that the DE needs.
     * return f1, double*/

    private static double getF1(DifferentialEquation de, double h, double x, double y){
        return h*de.evaluate(x,y);
    }

    /** getF2. this method returns the f2 value needed for Runge-Kutta methods
     * @param de - differential equation
     * @param h - growing interval of x
     * @param f1 - double f1 needed for f2 formula
     * @param x - x that the DE uses
     * @param y - y initial condition that the DE needs.
     * return f2, double*/

    private static double getF2(DifferentialEquation de, double h, double f1, double x, double y){
        return h*(de.evaluate(x+(1/2)*h,y+(1/2)*f1));
    }

    /** getF3. this method returns the f3 value needed for Runge-Kutta 4th order method.
     * @param de - differential equation
     * @param h - growing interval of x
     * @param f2 - double f2 needed for f3 formula
     * @param x - x that the DE uses
     * @param y - y initial condition that the DE needs.
     * return f3, double*/

    private static double getF3(DifferentialEquation de, double h, double f2, double x, double y){
        return h*(de.evaluate(x+(1/2)*h,y+(1/2)*f2));
    }

    /** getF4. this method returns the f4 value needed for Runge-Kutta 4th order method
     * @param de - differential equation
     * @param h - growing interval of x
     * @param f3 - double f3 needed for f4 formula
     * @param x - x that the DE uses
     * @param y - y initial condition that the DE needs.
     * return f4, double*/

    private static double getF4(DifferentialEquation de, double h, double f3, double x, double y){
        return h*(de.evaluate(x+h,y+f3));
    }

    /**Method that calculates the max of the second derivative in a given interval.
    * Needed for calculating the error in trapezium and Romberg integration methods.
    * @param f - Function to be used, for finding aproximations to its roots.
    * @param from - double that represents the first number to be evaluated, the initial point.
    * @param to - double that represents the end of the interval.
    * return max.
    * */

    private static double maxOfSecondDerivative(Function f, double from, double to){
        double max = f.secondDerivative(from);
        for(double i = from; i<=to;i+=0.001 ){
            double result = f.secondDerivative(i);
            if(result > max){
                max = result;
            }
        }
        return max;
    }

    /**Method that calculates the max of the fourth derivative in a given interval.
     * Needed for calculating the error in Simpson's integration method.
     * @param f - Function to be used, for finding aproximations to its roots.
     * @param from - double that represents the first number to be evaluated, the initial point.
     * @param to - double that represents the end of the interval.
     * return max.
     * */

    private static double maxOfFourthDerivative(Function f, double from, double to){
        double max = f.fourthDerivative(from);
        for(double i = from;i<=to;i+=0.001){
            double result = f.fourthDerivative(i);
            if(result>max){
                max = result;
            }
        }
        return max;
    }

    /**Method that checks if the sign of a number is positive or negative
     * @param a, number to be checked}
     * returns 1 if positive or zero, -1 if negative*/

    private static int sign(double a){
        if(a >= 0){
            return 1;
        }else{
            return -1;
        }
    }

    /**executeNewtonRaphson. Method that owns the logic of the newton raphson method described above.
     * @param f - Function f that newton raphson is using
     * @param x0 - initial value
     * return root*/

    private static double executeNewtonRaphson(Function f, double x0) {
        double xi = 0;
        try {
            xi = x0 - (f.evaluate(x0) / f.derivative(x0));
            if (xi - x0 != 0) {
               xi = executeNewtonRaphson(f, xi);
            }
            return xi;
        } catch (StackOverflowError s) {
            return xi;
        }
    }


    /**checkNRconditions. Method that checks if the conditions requested for the Newton-Raphson method are satisfied by
    * the function that is passed. The condition is:
    * f(from)*f(to)<0
    * f'(x) != 0 for every x belonging to [from,to]
    * f has no inflexion points.
    * @param f - Function that is going to be used,
    * @param from,to - limits of the interval [from,to]
    * return boolean
    * */

    private static boolean checkNRConditions(Function f, double from, double to){
        if(f.evaluate(from)*f.evaluate(to)<0){
            if(isNotZero(f,from,to)){
                if(hasNoInflexionPoint(f,from,to)){
                    return true;
                }
            }
        }else{
            return false;
        }
        return false;
    }

    /**isNotZero. condition needed for Newton-Raphson*/

    private static boolean isNotZero(Function f,double from, double to){
        for(double i = from;i<=to;i+=0.01){
            double result = f.derivative(i);
            if(result == 0){
                return false;
            }
        }
        return true;
    }

    /**hasNoInflexionPoint. condition needed for Newton-Raphson*/

    private static boolean hasNoInflexionPoint(Function f, double from, double to){
        ArrayList positive = new ArrayList();
        ArrayList negative = new ArrayList();
        int count1 = 0;
        for(double i = from;i<=to;i+=0.2){
            double result = f.secondDerivative(i);
            if(result > 0){
                positive.add(result);
            }else if(result < 0){
                negative.add(result);
            }else{
                if(positive.size()>negative.size()){
                    positive.add(result);
                }else{
                    negative.add(result);
                }
            }
            count1++;
        }
        if(count1 == positive.size()){
            return true;
        }else if(count1 == negative.size()){
            return true;
        }else{
            return false;
        }
    }

}




