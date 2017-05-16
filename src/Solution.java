public class Solution {
    private static double calcF1(double x, double y) {
        return Math.pow(x, 2) / 4 + Math.pow(y, 2) / 9 - 1;
    }

    private static double calcF2(double x, double y) {
        return x - Math.pow(y, 2);
    }

    private static double calcF1DerivX(double x, double y) {
        return x / 2;
    }

    private static double calcF1DerivY(double x, double y) {
        return 2 * y / 9;
    }

    private static double calcF2DerivX(double x, double y) {
        return 1;
    }

    private static double calcF2DerivY(double x, double y) {
        return -2 * y;
    }

    private static double calcFi1DerivX(double x, double y, double h) {
        return (calcF1(x, y) - calcF1(x - h, y)) / h;
    }

    private static double calcFi1DerivY(double x, double y, double h) {
        return (calcF1(x, y) - calcF1(x, y - h)) / h;
    }

    private static double calcFi2DerivX(double x, double y, double h) {
        return (calcF2(x, y) - calcF2(x - h, y)) / h;
    }

    private static double calcFi2DerivY(double x, double y, double h) {
        return (calcF2(x, y) - calcF2(x, y - h)) / h;
    }

    private static double epsilon = 0.00001;

    public static void main(String[] args) {
        double xcur;
        double xnext = 1.3026489924004072;
        double ycur;
        double ynext = 1;
        double h1 = 0.0000001;
        double h2 = 0.0000001;
        double a;
        int iterations = 0;
        System.out.println("X приближённое: " + xnext);
        System.out.println("Y приближённое: " + ynext);
        do {
            xcur = xnext;
            ycur = ynext;
            a = 1 / (calcF1DerivX(xcur, ycur) * calcF2DerivY(xcur, ycur) - calcF1DerivY(xcur, ycur) * calcF2DerivX(xcur, ycur));
            xnext = xcur - a * (calcFi2DerivY(xcur, ycur, h2) * calcF1(xcur, ycur) - calcFi1DerivY(xcur, ycur, h2) * calcF2(xcur, ycur));
            ynext = ycur - a * (calcFi1DerivX(xcur, ycur, h1) * calcF2(xcur, ycur) - calcFi2DerivX(xcur, ycur, h1) * calcF1(xcur, ycur));
            iterations++;
        } while (Math.max(Math.abs(xnext - xcur), Math.abs(ynext - ycur)) >= epsilon);
        System.out.println("X найденное: " + xnext);
        System.out.println("Y найденное: " + ynext);
        System.out.println("Невязка f1: " + Math.abs(calcF1(xnext, ynext)));
        System.out.println("Невязка f2: " + Math.abs(calcF2(xnext, ynext)));
        System.out.println("Количество итераций: " + iterations);
    }
}
