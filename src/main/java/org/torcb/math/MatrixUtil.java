package org.torcb.math;

import static org.torcb.math.SymbolMath.*;

public interface MatrixUtil {

    @SuppressWarnings("unused")
    class AxisToRot3D {
        double i, j, k; //|ijk|=1
        double angleRadians;
        Matrix L;
        double sin, cos;
        Matrix LL, rot3d;
        double[] rot3dScalar;

        AxisToRot3D(double u, double v, double w) {
            this(null, u, v, w);
        }

        AxisToRot3D(Double rotAngR, double u, double v, double w) {
            var lenUnitVec = LenUnitVector.from(u, v, w);
            if (zero(lenUnitVec.len)) return;
            double len = lenUnitVec.len;
            i = lenUnitVec.unitVector[0];
            j = lenUnitVec.unitVector[1];
            k = lenUnitVec.unitVector[2];
            L = new Matrix(3, 3).init(
                    0, -k, j,
                    k, 0, -i,
                    -j, i, 0);
            LL = L.multiplyIm(L);
            angleRadians = rotAngR == null ? len : rotAngR;
            sin = Math.sin(angleRadians);
            cos = Math.cos(angleRadians);
            rot3d = Matrix.identity(3)
                          .addIm(L.multiplyIm(sin))
                          .addIm(LL.multiplyIm(1.0 - cos));
            rot3dScalar = rot3d.getAllScalars();
            //rot3d      = R = I + L*sin + L*L*(1-cos)
            //invRot_ijk = L = (R - tr(R)) / (2*sin)
        }

        public Matrix getRotateMatrix() {
            return rot3d;
        }

        public double[] getRotateScalar() {
            return rot3dScalar;
        }
    }

    static double length(double[] vector) {
        return Math.sqrt(dotProduct(vector, vector));
    }

    record LenUnitVector(double len, double[] unitVector) {
        public static LenUnitVector from(double[] vector, int... indices) {
            double[] slice = new double[indices.length];
            int i = 0;
            for (int ix : indices) {
                slice[i++] = vector[ix];
            }
            return from(slice);
        }

        public static LenUnitVector from(double... vector) {
            var out = new LenUnitVector(length(vector), new double[vector.length]);
            if (out.len > 0) {
                for (int i = 0; i < vector.length; i++) {
                    out.unitVector[i] = vector[i] / out.len;
                }
            }
            return out;
        }
    }

    static double dotProduct(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            sum += a[i] * b[i];
        }
        return sum;
    }

    @SuppressWarnings("unused")
    class Rot3dToAxis {
        final double[] rot;
        final double traceSum;
        final Matrix rotRotTr;
        final double[] rotRotTrSc;
        final LenUnitVector lenUnitVector;
        final double cos_fi, sin_fi;
        final double angleRadians;

        public Rot3dToAxis(double[] rotate3d) {
            if (rotate3d.length != 9)
                throw new IllegalArgumentException("Expects 3x3 values");
            this.rot = rotate3d;
            traceSum = rotate3d[0] + rotate3d[4] + rotate3d[8];
            rotRotTr = new Matrix(3, 3).init(rotate3d).op(
                    r -> r.minusIm(r.transposeIm())
            );
            rotRotTrSc = rotRotTr.getAllScalars();
            lenUnitVector = LenUnitVector.from(rotRotTrSc, 7, 2, 3);
            cos_fi = (traceSum - 1.0) / 2.0; //cos = (trace-1)/2
            sin_fi = lenUnitVector.len / 2.0;
            angleRadians = Math.atan2(sin_fi, cos_fi);
        }

        public double[] getRotationAxis() {
            return lenUnitVector.unitVector;
        }

        public double getRotationAngleRadians() {
            return angleRadians;
        }
    }

    record Quaternion(MVPolynomial c, MVPolynomial i, MVPolynomial j, MVPolynomial k) {
        public static Quaternion fromSymbols(String c, String i, String j, String k) {
            return new Quaternion(
                    MVPolynomial.parse(c),
                    MVPolynomial.parse(i),
                    MVPolynomial.parse(j),
                    MVPolynomial.parse(k)
            );
        }

        public MVPolynomial[] asArr() {return new MVPolynomial[]{c, i, j, k};}

        public Quaternion conjugate() {return new Quaternion(c, i.negateIm(), j.negateIm(), k.negateIm());}

        public Quaternion multiply(Quaternion o) {
            var _i = i.negateIm();
            var _j = j.negateIm();
            var _k = k.negateIm();
            MVPolynomial[] other = o.asArr();
            return new Quaternion(
                    SymbolMath.dotProduct(other, c, _i, _j, _k),
                    SymbolMath.dotProduct(other, i, c, k, _j),
                    SymbolMath.dotProduct(other, j, _k, c, i),
                    SymbolMath.dotProduct(other, k, j, _i, c)
            );
        }

        public Quaternion substituteTerms(Term term, MVPolynomial toExpr) {
            return new Quaternion(
                    c.substituteTermsIm(term, toExpr),
                    i.substituteTermsIm(term, toExpr),
                    j.substituteTermsIm(term, toExpr),
                    k.substituteTermsIm(term, toExpr)
            );
        }

        public Quaternion substituteTerms(String s, String expr) {
            var term = new Term(s);
            var toExpr = MVPolynomial.parse(expr);
            return substituteTerms(term, toExpr);
        }

        @Override
        public String toString() {
            return "Quaternion:\n c:  " + c +
                   "\n i:  " + i +
                   "\n j:  " + j +
                   "\n k:  " + k;
        }
    }
}
