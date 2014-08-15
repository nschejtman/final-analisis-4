package function.knownFunction;

import function.Function;

public class LinearFunction implements Function {
    double a, b;

    public LinearFunction(double a, double b) {
        this.a = a;
        this.b = b;
    }

    @Override
    public double evaluate(double x){
        return (a*x)+b;
    }

    public double root(){
        return -(b/a);
    }

    @Override
    public Function derivative() {
        return new ConstantFunction(a);
    }

    @Override
    public Function integral() {
        return null;  //TODO
    }
}
