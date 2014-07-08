package finalManu;

/**
 * Created by Manuel on 23/02/14.
 */
public class Cuadratic extends Function {

    private double a,b,c;

    public Cuadratic(double a,double b, double c){
        this.a = a;
        this.b = b;
        this.c = c;
    }

    @Override
    double evaluate(double x) {
        return a*Math.pow(x,2) + b*x + c;
    }

    @Override
    double derivative(double x) {
        return 2*a*x + b;
    }

    @Override
    double secondDerivative(double x) {
        return 2*a;
    }

    @Override
    double thirdDerivative(double x) {
        return 0;
    }

    @Override
    double fourthDerivative(double x) {
        return 0;
    }
}
