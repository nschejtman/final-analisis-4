package function.knownFunction;

import function.Function;

public class ConstantFunction implements Function {
    double value;

    public ConstantFunction(double value) {
        this.value = value;
    }

    public double getValue() {
        return value;
    }

    public void setValue(double value) {
        this.value = value;
    }

    @Override
    public Function derivative() {
        return null;
    }

    @Override
    public Function integral() {
        return new LinearFunction(value, Math.random());
    }

    @Override
    public double evaluate(double x) {
        return value;
    }
}
