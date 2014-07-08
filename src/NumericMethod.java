public abstract class NumericMethod {

    /**
     * Bisection numeric method to approximate the root of a f(x) function in a given interval. This method has three
     * conditions:
     * 1. There should be only one root in the interval
     * 2. f(from) * f(to) < 0 -> has a root in the interval
     * 3. df/dx != 0 -> is not constant
     *
     * @param function function to be analyzed
     * @param from     lower bound
     * @param to       upper bound
     * @return root
     */
    public static double bisection(Function function, double from, double to) {
        //check df/dx != 0
        Function derivative = function.derivative();
        boolean isConstant = true;
        for (int i = 0; i < 500; i++) {
            if (derivative.evaluate(i) != 0) isConstant = false;
        }
        if (isConstant) throw new IllegalArgumentException("Function is constant");

        //check f(from)*f(to) != 0
        final double yFrom = function.evaluate(from);
        final double yTo = function.evaluate(to);
        if(yFrom*yTo >= 0) throw new IllegalArgumentException("None or multiple roots in the interval");

        //Apply method
        final double middle = (from + to)/2;
        final double yMiddle = function.evaluate(middle);

        if(yMiddle == 0) return middle;
        else{
            if(yFrom*yMiddle < 0) return bisection(function, from, middle);
            else return bisection(function, middle, to);
        }
    }

    public static double newtonRaphson(Function function, double x0){
        if(function.evaluate(x0) == 0) return 0;
        else{
            double x1 = x0 - (function.evaluate(x0)/function.derivative().evaluate(x0));
            return newtonRaphson(function, x1);
        }
    }


}
