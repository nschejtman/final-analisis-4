import function.Function;
import function.knownFunction.FunctionTest;
import function.knownFunction.LinearFunction;
import function.knownFunction.PolynomialFunction;
import function.knownFunction.QuadraticFunction;
import function.unknownFunction.DiffEqImpl;

public class Main {
  public static void main(String args[])
  {
//     QuadraticFunction quadraticFunction = new QuadraticFunction(2,2,2);
//
////
//      System.out.println("SIMPSON" + NumericMethod.simpson(new FunctionTest(),0,3.14,0.01));
//      System.out.println("Trapecio" + NumericMethod.trapeziums(new FunctionTest(), 0, 3.14, 0.01));
//      System.out.println("Romberg"+ NumericMethod.integrateByRomberg(new FunctionTest(),0,3.14,2,0));
////
//      LinearFunction linearFunction = new LinearFunction(1, -1);
//      PolynomialFunction polynomialFunction = new PolynomialFunction(2);
//      polynomialFunction.setCoefficient(2,2);
//      polynomialFunction.setCoefficient(1,2);
//      polynomialFunction.setCoefficient(0,2);
//      System.out.println(NumericMethod.newtonRaphson(linearFunction, 0, 2));
//      System.out.println("Trapecio2: " + NumericMethod.trapeziums(polynomialFunction, 1, 3, 0.001));

      //Function function = new QuadraticFunction(1, -3 , 2);
     // System.out.println(NumericMethod.newtonRaphson(function, 0.5, 1.6, 0.001));

      DiffEqImpl dq = new DiffEqImpl();

      System.out.println(NumericMethod.rungeKuttaO4(dq, 0, 1, 0.1, 2));


  }
}
