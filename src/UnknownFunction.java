/**
 * Class for an unknown function. As the function is unknown, its derivative and integral are as well.
 */
public abstract class UnknownFunction implements Function {
    @Override
    public Function derivative() {
        return null;
    }

    @Override
    public Function integral() {
        return null;
    }

    @Override
    public abstract double evaluate(double x);

    public abstract double firstDerivative(double x);
    public abstract double secondDerivative(double x);
    public abstract double thirdDerivative(double x);
    public abstract double fourthDerivative(double x);
}
