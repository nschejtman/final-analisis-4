package function.knownFunction;

import function.Function;

public class FunctionTest2 implements Function {
    @Override
    public Function derivative() {
        return null;
    }

    @Override
    public Function integral() {
        return null;
    }

    @Override
    public double evaluate(double x) {
        return 1 - Math.sin(x);
    }
}
