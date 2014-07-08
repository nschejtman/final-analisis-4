package finalManu;

/**
 * Created by Manuel on 18/02/14.
 */
public class Main {

    public static void main(String[] args) {
        Function f = new Function() {

            @Override
            double evaluate(double x) {
                double result = Math.pow(x,3) - 2*Math.pow(x,2) -5;
                return result;
            }

            @Override
            double derivative(double x) {
                double result = 3*Math.pow(x,2) -4*x;
                return result;
            }

            @Override
            double secondDerivative(double x) {
                double result = 6*x - 4;
                return result;
            }

            @Override
            double thirdDerivative(double x) {
                return 0;
            }

            @Override
            double fourthDerivative(double x) {
                return 0;
            }
        };
        Function f2 = new Function() {

            @Override
            double evaluate(double x) {
                double result = x - Math.cos(x);
                return result;
            }

            @Override
            double derivative(double x) {
                double result = 1 + Math.sin(x);
                return result;
            }

            @Override
            double secondDerivative(double x) {
                double result = Math.cos(x);
                return result;
            }
            @Override
            double thirdDerivative(double x) {
                return 0;
            }

            @Override
            double fourthDerivative(double x) {
                return 0;
            }

        };
        Function f3 = new Function() {

            @Override
            double evaluate(double x) {
               return (Math.pow(x,2))-2;
            }

            @Override
            double derivative(double x) {
                return 2*x;
            }

            @Override
            double secondDerivative(double x) {
                return 0;
            }

            @Override
            double thirdDerivative(double x) {
                return 0;
            }

            @Override
            double fourthDerivative(double x) {
                return 0;
            }
        };
        Function f4 = new Function() {

            @Override
            double evaluate(double x) {
                return 1/(x);
            }

            @Override
            double derivative(double x) {
                return -1/(Math.pow(x,2));
            }

            @Override
            double secondDerivative(double x) {
                return  2/(Math.pow(x+1,3));
            }

            @Override
            double thirdDerivative(double x) {
                return 0;
            }

            @Override
            double fourthDerivative(double x) {
                return 0;
            }
        };
        Function f5 = new Function() {
            @Override
            double evaluate(double x) {
                return (Math.sin(x))/x;
            }

            @Override
            double derivative(double x) {
                return 0;
            }

            @Override
            double secondDerivative(double x) {
                return 0;
            }

            @Override
            double thirdDerivative(double x) {
                return 0;
            }

            @Override
            double fourthDerivative(double x) {
                return 0;
            }
        };
        Cuadratic c = new Cuadratic(2,4,5);
        DifferentialEquation de = new DifferentialEquation() {
            @Override
            double evaluate(double x, double y) {
                return Math.cos(x) - Math.sin(y) + Math.pow(x,2);
            }

            @Override
            double derivative(double x, double y) {
                return -Math.sin(x) - Math.cos(y)*evaluate(x,y) + 2*x;
            }

            @Override
            double secondDerivative(double x, double y) {
                return -Math.cos(x) + Math.sin(y)*Math.pow(derivative(x,y),2) - Math.cos(y)*derivative(x, y) + 2;
            }

            @Override
            double thirdDerivative(double x, double y) {
                return 0;
            }

            @Override
            double fourthDerivative(double x, double y) {
                return 0;
            }
        };

       //Numeric.newtonRaphson(f2, 0, (Math.PI)/2);
       //Numeric.bisection(f3,0,2);
       Numeric.integrateByTrapezium(f,1,3,0.1);
        System.out.println(Numeric.integrateBySimpson(f,1,3,0.1));
       //Numeric.integrateByRomberg(f,1,3,2,0);
       //Numeric.solveDiffEqByTaylor(de,-1,1,0.01,3);
       //Numeric.rungeKutta4thOrder(de,-1,1,0.01,3);


    }

}
