/**
 * Linear function of type f(x) = a x + b
 */
public class LinearFunction {
    private double a;
    private double b;

    public LinearFunction(double a, double b) {
        this.a = a;
        this.b = b;
    }

    public double evaluate(double x){
        return (a*x)+b;
    }

    public double root(){
        return -(b/a);
    }


}
