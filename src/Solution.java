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

    private static double calcMatrixNorm(double[][] mtr) {
        double result = -1;
        double sum;
        for(int i = 0; i < mtr.length; i++) {
            sum = 0;
            for (int j = 0; j < mtr[i].length; j++) {
                sum += Math.abs(mtr[i][j]);
            }
            result = Math.max(result, sum);
        }
        return result;
    }

    private static double calcVectorNorm(double[] vector) {
        double result = Math.abs(vector[0]);
        for(int i = 1; i < vector.length; i++) {
            result = Math.max(result, vector[i]);
        }
        return result;
    }

    private static double[] calcVectorSubtract(double[] left, double[] right) {
        double[] result = left.clone();
        for (int i = 0; i < left.length; i++) {
            result[i] -= right[i];
        }
        return result;
    }

    private static double[] calcMatrixMultiplyVector(double[][] mtr, double[] vector) {
        double[] result = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            result[i] = 0;
            for (int j = 0; j < mtr[i].length; j++) {
                result[i] += mtr[i][j] * vector[j];
            }
        }
        return result;
    }

    private static double[][] calcJacobian(double x, double y) {
        double[][] result = new double[2][2];
        result[0][0] = calcF1DerivX(x, y);
        result[0][1] = calcF1DerivY(x, y);
        result[1][0] = calcF2DerivX(x, y);
        result[1][1] = calcF2DerivY(x, y);
        return result;
    }

    private static double[][] calcInverseJacobian(double x, double y) {
        double[][] result = new double[2][2];
        double c = 1 / (calcF1DerivX(x, y) * calcF2DerivY(x, y) - calcF1DerivY(x, y) * calcF2DerivX(x, y));
        result[0][0] = c * calcF2DerivY(x, y);
        result[0][1] = -c * calcF1DerivY(x, y);
        result[1][0] = -c * calcF2DerivX(x, y);
        result[1][1] = c * calcF1DerivX(x, y);
        return result;
    }

    private static void checkTheoremConditions(double x, double y) {
        double[][] jacobian;
        double[][] jacobianInverse;
        double[] fVect = new double[2];
        double[] exactSolution = new double[2];
        double[] approxSolution = new double[2];
        double a1;
        double a2;
        double delta;
        exactSolution[0] = 1.7900855862527592;
        exactSolution[1] = 1.3379408007280291;
        approxSolution[0] = x;
        approxSolution[1] = y;
        delta = calcVectorNorm(calcVectorSubtract(exactSolution, approxSolution));
        fVect[0] = calcF1(exactSolution[0], exactSolution[1]) - calcF1(x, y);
        fVect[1] = calcF2(exactSolution[0], exactSolution[1]) - calcF2(x, y);
        jacobian = calcJacobian(exactSolution[0], exactSolution[1]);
        jacobianInverse = calcInverseJacobian(exactSolution[0], exactSolution[1]);
        a1 = calcMatrixNorm(jacobianInverse);
        a2 = calcVectorNorm(calcVectorSubtract(fVect, calcMatrixMultiplyVector(jacobian, calcVectorSubtract(exactSolution, approxSolution)))) /
                Math.pow(calcVectorNorm(calcVectorSubtract(exactSolution, approxSolution)), 2);
        System.out.println("a1: " + a1);
        System.out.println("a2: " + a2);
        System.out.println("c: " + a1 * a2);
        System.out.println("delta: " + delta);
        System.out.println("b: " + Math.min(1 / (a1 * a2), delta));
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
        checkTheoremConditions(xnext, ynext);
        System.out.println("X приближённое: " + xnext);
        System.out.println("Y приближённое: " + ynext);
        do {
            xcur = xnext;
            ycur = ynext;
            a = 1 / (calcFi1DerivX(xcur, ycur, h1) * calcFi2DerivY(xcur, ycur, h2) - calcFi1DerivY(xcur, ycur, h1) * calcFi2DerivX(xcur, ycur, h2));
            xnext = xcur - a * (calcFi2DerivY(xcur, ycur, h2) * calcF1(xcur, ycur) - calcFi1DerivY(xcur, ycur, h2) * calcF2(xcur, ycur));
            ynext = ycur - a * (calcFi1DerivX(xcur, ycur, h1) * calcF2(xcur, ycur) - calcFi2DerivX(xcur, ycur, h1) * calcF1(xcur, ycur));
            iterations++;
        } while (Math.max(Math.abs(xnext - xcur), Math.abs(ynext - ycur)) >= epsilon);
        System.out.println("X найденное: " + xnext);
        System.out.println("Y найденное: " + ynext);
        System.out.println("Невязка: " + Math.max(Math.abs(calcF1(xnext, ynext)), Math.abs(calcF2(xnext, ynext))));
        System.out.println("Количество итераций: " + iterations);
    }
}
