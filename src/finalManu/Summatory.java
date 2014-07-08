package finalManu;

/**
 * Created by Manuel on 19/02/14.
 */
public class Summatory {

    /*Class Summatory. This class has static methods that return the value of different summatories,
      results that are needed for integration by romberg, simpson and trapecios.
      @author: Manuel Acuña
      @author: Nicolás Rudolph
    * */
    private Summatory(){};

    /*trapecios: returns the value of the summatory used in trapecios's formula. the summatory is:
      the summatory of f(xi) from i = 1 to i = n-1, being xi = x0 + i*h, and n = (xn - x0)/h. h refers to
      the growing interval.
      @param Function f - the function that the user desires to integrate
      @param double a - initial value, value that is going to be used by the function, in order to evaluate itself and its derivatives
      @param double h - growing interval.
      @param double n - (b-a)/h
     *  */

    public static double trapecios(Function f, double a, double h, double n){
        double summatory = 0;
        double xi;
        for(int i = 1;i<n;i++){
            xi = a+i*h;
            summatory = summatory + f.evaluate(xi);
        }
        return summatory;
    }

    /*simpsonOdd: returns the value of the first summatory in Simpson formula for integration. The summatory is:
      Summatory of(f(x2j)) from j = 1 to j = (n/2)-1; being x2j = x0 + 2*j*h, and n = (xn-x0)/h. h refers to
      the growing interval.
      @param Function f - the function that the user desires to integrate
      @param double a - initial value, value that is going to be used by the function, in order to evaluate itself and its derivatives
      @param double h - growing interval.
      @param double n - (b-a)/h
     *  */

    public static double simpsonOdd(Function f, double a, double h, double n){
        double summatory = 0;
        double x2j;
        for(int j = 1;j<=((n/2)-1);j++){
            x2j = a + 2*j*h;
            summatory = summatory + f.evaluate(x2j);
        }
        System.out.println("odd: " +summatory);
        return summatory;
    }

    /*simpsonEven: returns the value of the second summatory inSimpson formula for integration. The summatory is:
    * Summatory of(f(x2j - 1)) from j = 1 to (n/2); being x2j-1 = x0 + (2*j-1)*h, and n = (xn-x0)/h. h refers to
      the growing interval.
      @param Function f - the function that the user desires to integrate
      @param double a - initial value, value that is going to be used by the function, in order to evaluate itself and its derivatives
      @param double h - growing interval.
      @param double n - (b-a)/h }
      */

    public static double simpsonEven(Function f, double a, double h, double n){
        double summatory = 0;
        double xEven;
        for(int j = 1;j<=(n/2);j++){
            xEven = a + (2*j-1)*h;
            summatory = summatory + f.evaluate(xEven);
        }
        System.out.println("even: "+summatory);
        return summatory;
    }

}
