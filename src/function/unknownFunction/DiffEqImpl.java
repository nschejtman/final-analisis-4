package function.unknownFunction;

public class DiffEqImpl extends DifferentialEquation {
    @Override
    public double evaluate(double x, double y) {
        return 3*x*y;
    }

    @Override
    public double derivative(double x, double y) {
        return 3*y + 3*x*evaluate(x, y);
    }

    @Override
    public double secondDerivative(double x, double y) {
        return 6;
    }

    @Override
    public double thirdDerivative(double x, double y) {
        return 0;
    }

    @Override
    public double fourthDerivative(double x, double y) {
        return 0;
    }
}
