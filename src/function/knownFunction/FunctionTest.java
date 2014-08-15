package function.knownFunction;

import function.Function;

public class FunctionTest implements Function {


    @Override
    public Function derivative() {
        return new FunctionTest2();
    }

    @Override
    public Function integral() {
        return null;
    }

    @Override
    public double evaluate(double x){
        return x - Math.cos(x);
    }
}
