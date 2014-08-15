package function.knownFunction;

import function.Function;

public class QuadraticFunction implements Function {
    double a, b, c;

    public QuadraticFunction(double a, double b, double c) {
        this.a = a;
        this.b = b;
        this.c = c;
    }

    public double getA() {
        return a;
    }

    public double getB() {
        return b;
    }

    public double getC() {
        return c;
    }

    @Override
    public Function derivative() {

        return new LinearFunction(a, b);
    }

    @Override
    public Function integral() {
        PolynomialFunction integ = new PolynomialFunction(3);
        integ.setCoefficient(0, Math.random());
        integ.setCoefficient(1, c);
        integ.setCoefficient(2, b);
        integ.setCoefficient(3, a);
        return integ;
    }

    @Override
    public double evaluate(double x) {
        return a * Math.pow(x, 2) + b * x + c;
    }
}
