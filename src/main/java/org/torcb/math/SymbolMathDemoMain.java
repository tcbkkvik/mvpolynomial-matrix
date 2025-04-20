package org.torcb.math;

import java.util.Arrays;
import java.util.function.BiFunction;

import static org.torcb.math.SymbolMath.*;
import static org.torcb.math.MatrixUtil.Quaternion;

@SuppressWarnings("CommentedOutCode")
public interface SymbolMathDemoMain {
    static void tstEulerMatrixPow() {
        System.out.println("------------------------------");
        System.out.println(" tstEulerMatrixPow");
        System.out.println("------------------------------");
        var l90_1 = Matrix.init3x3(
                "0", "-k", "j",
                "k", "0", "-i",
                "-j", "i", "0").label("xProd3x3");

        var iiRule = SubstituteTerm.parse("i i", "1 - j j - k k");
        var jjRule = SubstituteTerm.parse("j j", "1 - i i - k k");
        var l90_2 = l90_1.multiplyIm(l90_1)
                         .substituteTermsIm(jjRule)
                         .substituteTermsIm(iiRule)
                         .substituteTermsIm(jjRule);
        var l90_3 = l90_2.multiplyIm(l90_1)
                         .substituteTermsIm(iiRule);
        var l90_4 = l90_2.multiplyIm(l90_2)
                         .substituteTermsIm(iiRule)
                         .substituteTermsIm(jjRule);
        Matrix neg_3 = l90_3.multiplyIm(-1);

        var r1_3 = l90_1.addIm(l90_3);
        var r2_4 = l90_2.addIm(l90_4);

        var zeroMatrix = Matrix.init3x3(
                "0", "0", "0",
                "0", "0", "0",
                "0", "0", "0");
        boolean isNull_1 = r1_3.equals(zeroMatrix);
        boolean isNull_2 = r2_4.equals(zeroMatrix);

        System.out.println("l90_1" + l90_1);
        System.out.println("l90_2" + l90_2);
        System.out.println("l90_3" + l90_3);
        System.out.println("l90_4" + l90_4);
        System.out.println("add 1 3" + r1_3);
        System.out.println("add 2 4" + r2_4);
        System.out.println("add 1 3 ==0 " + isNull_1);
        System.out.println("add 2 4 ==0 " + isNull_2);
        assert l90_1.equals(neg_3);
        assert isNull_1;
        assert isNull_2;
    }

    static void tstEulerRotMatrixTransposeProdIsIdentity() {
        System.out.println("------------------------------");
        System.out.println(" tstEulerRotMatrixTransposeProdIsIdentity");
        System.out.println("------------------------------");
        /* https://thenumb.at/Exponential-Rotations/#euler-angles
            rotMat = L*L − L*L*cos + I + L*sin
            transpose(rotMat) * rotMat = I
        */
        //product of 3x3 rotation, transpose(R)*R
        var rotMat = MVPolynomial
                .parse("1 - cos")
                .multiplyIm("L L")
                .add("I + L sin");
        var rotMatTranspose = MVPolynomial
                .parse("1 - cos")
                .multiplyIm(" L L")
                .add("I - L sin");
        var prod_1 = rotMatTranspose.multiplyIm(rotMat);
        var prod_2 = prod_1.substituteTermsIm("L L L L", "- L L");
        var prod_3 = prod_2.substituteTermsIm("sin sin", "1 - cos cos");
        var prod_4 = prod_3.substituteTermsIm("I L", "L");
        var prod_5 = prod_4.substituteTermsIm("I I", "I");

        Matrix.printMatrixRingBufAndClear();

        //Result = Identity matrix 'I' expected
        System.out.println("tr(rotMat)= " + rotMatTranspose);
        System.out.println("rotMat = " + rotMat);
        System.out.println("tr(rotMat) * rotMat");
        System.out.println(" = " + prod_1);
        System.out.println(" = " + prod_2);
        System.out.println(" = " + prod_3);
        System.out.println(" = " + prod_4);
        System.out.println(" = " + prod_5);
    }

    private static void testMisc() {
        System.out.println("------------------------------");
        System.out.println(" testMisc");
        System.out.println("------------------------------");
        var m21 = new Matrix(2, 1).init("a", "b").label("m21");
        var m12 = m21.transposeIm().label("m12");
        Matrix.RowCol rowCol = (pos, row, col, val) ->
                System.out.printf(" pos=%d, row=%d, col=%d, val=%s\n", pos, row, col, val);
        m21.iterate(rowCol);
        System.out.println();
        m12.iterate(rowCol);
//        System.out.println("m21:" + m21);
//        System.out.println("tr(m21):" + m12);

        var halfRotate = Quaternion.fromSymbols("cos", "sin", null, null);
        var halfRotateCnj = halfRotate.conjugate();
        var qVector = Quaternion.fromSymbols(null, "x", "y", "z");
        var v1 = halfRotate.multiply(qVector);
        var v2 = v1.multiply(halfRotateCnj);
        System.out.println(qVector);
        System.out.println(v1);
        System.out.println(v2);

        var v3a = v2.substituteTerms("sin sin", "ss")
                    .substituteTerms("cos sin", "cs")
                    .substituteTerms("cos cos", "cc");
        System.out.println("v3a");
        System.out.println(v3a);

        // Rule:  sin*sin = 1 - cos*cos
        var v3 = v2.substituteTerms("sin sin", "1 - cos cos");
        var v4 = v3.substituteTerms("sin cos", "0.5 sin2")
                   .substituteTerms("cos cos", "0.5 + 0.5 cos2");
        System.out.println("After substituteTerms -> v3:");
        System.out.println(v3);
        System.out.println(v4);
        System.out.println();

        var a = new MVPolynomial()
                .add(1, "sin")
                .add(2, "x")
                .add(1, "y");
        var b = new MVPolynomial()
                .add(2, "x")
                .add(-1, "y");
        var c = a.multiplyIm(b);
        System.out.println(" a:        " + a);
        System.out.println(" b:        " + b);
        System.out.println(" a*b=c:    " + c);
    }

    @SuppressWarnings("unused")
    static Matrix calcRot3x3Matrix() {
        Matrix.logRingBuf.clear();
        SUBSTITUTE_RULES.set(new SubstituteTerms()
                .add("i i", "1 - j j - k k")
                .add("cos cos", "1 - sin sin"));
        var L = Matrix.init3x3(
                "0", "-k", "j",
                "k", "0", "-i",
                "-j", "i", "0").label("L (CrossProduct matrix; Lie algebra)");
        var L_sin = L.multiplyIm("sin").label("L*sin");
        var LL = L.multiplyIm(L).label("L*L");
        var LL_1_cos = LL.multiplyIm("1 - cos").label("LL*(1-cos)");// 1 - cos;
        Matrix rotateM = Matrix.identity(3).addIm(L_sin).addIm(LL_1_cos)
                               .label("Rotate3D = I + L*sin - L*L*(1-cos)  //Rodrigues formula");
        // - R*tr(T)==I:
        var idRot = rotateM.multiplyIm(rotateM.transposeIm()).label("I = Rotate3D * transpose(Rotate3D)");
        Matrix.printMatrixRingBufAndClear();

        var L_LTrans = L.multiplyIm(L.transposeIm());
        System.out.println("L*tr(L)" + L_LTrans);
        Matrix rotateM_deriveSin = rotateM.deriveIm("sin");
        var ll_derive_i = LL.deriveIm("i");
        Matrix rotateM_deriveCos = rotateM.deriveIm("cos");
        Matrix minusTrRotateM = rotateM.transposeIm().multiplyIm(-1);

        System.out.println("rotateM_deriveSin" + rotateM_deriveSin);
        System.out.println("rotateM_deriveCos" + rotateM_deriveCos);
        System.out.println("derive LL(i)" + ll_derive_i);
        System.out.println("neg_rotateM_tr" + minusTrRotateM);

        var rMinus_rTr = rotateM.addIm(minusTrRotateM);
        var invRot_ijk = rMinus_rTr.substituteTermsIm("sin", "0.5");

        MVPolynomial[] cells = invRot_ijk.cells;
        MVPolynomial[] ijk = {cells[7], cells[2], cells[3]};

        System.out.println("rMinus_rTr" + rMinus_rTr);
        System.out.println("invRot_ijk" + invRot_ijk);
        System.out.println("Restored original rotation_vector[ijk] = " + Arrays.toString(ijk));
        Matrix.logRingBuf.clear();
        SUBSTITUTE_RULES.remove();
        return rotateM;
    }

    static void runCalcRot3x3Matrix() {
        System.out.println("------------------------------");
        System.out.println(" runCalcRot3x3Matrix");
        System.out.println("------------------------------");

        double i = .7;
        double j = .2;
        double k = .1;
        var axis_toRot = new MatrixUtil.AxisToRot3D(i, j, k);
        var rot = axis_toRot.rot3d;

        Matrix L2 = rot.minusIm(rot.transposeIm())
                       .multiplyIm(1.0 / (2 * axis_toRot.sin));
        MVPolynomial[] ijk_vector = {
                L2.cell(2, 1),
                L2.cell(0, 2),
                L2.cell(1, 0)
        };
        var ijk_faster = new Matrix(1, 3).init(
                rot.cell(2, 1).minusIm(rot.cell(1, 2)),
                rot.cell(0, 2).minusIm(rot.cell(2, 0)),
                rot.cell(1, 0).minusIm(rot.cell(0, 1))
        ).multiplyIm(1.0 / (2 * axis_toRot.sin));

        MVPolynomial cosineAngleST = Matrix.trace(rot).add(-1).multiplyIm(0.5);
        double cosine_angle = cosineAngleST.scalarSum();
        double angleRad_recalculated = Math.acos(cosine_angle);
        double angle_err = angleRad_recalculated - axis_toRot.angleRadians;

        System.out.println("L orig:" + axis_toRot.L);
        System.out.println("L recreated:" + L2);
        System.out.println("L -> ijk_vector: " + Arrays.toString(ijk_vector));
        System.out.println("L -> ijk_faster: " + ijk_faster);
        System.out.println("cosineAngleST = MatrixUtil.trace(rot).add(-1).multiplyIm(0.5) = " + cosineAngleST);
        System.out.println("angleRad_recalculated = " + angleRad_recalculated);
        System.out.println("angle_err = " + angle_err);
        var rot_toAxis = new MatrixUtil.Rot3dToAxis(rot.getAllScalars());
        System.out.println("rot_toAxis.angleRadians: " + rot_toAxis.angleRadians);

        var rotTr = rot.transposeIm();
        var rotTrRot = rotTr.multiplyIm(rot);
        var identity = Matrix.identity(3);
        boolean idEqualOk = identity.equals(rotTrRot);

        System.out.printf("i,j,k: %f %f %f\n", i, j, k);
        System.out.println("-> rot" + rot);
        System.out.println("tr(rot)*rot" + rotTrRot);
        System.out.println("idEqualOk: " + idEqualOk);

        assert idEqualOk;
    }

    static void diff() {
        System.out.println("------------------------------");
        System.out.println(" derive, integrate");
        System.out.println("------------------------------");
        var poly = MVPolynomial.parse("3 x x a + 2 x b+ 9 c");
        var expect = MVPolynomial.parse("2b + 6a*x");
        var derived = poly.deriveIm("x");
        System.out.println("Polynomial f(x) = " + poly);
        System.out.println("Polynomial f`(x) = " + derived);
        boolean ok = expect.equals(derived);
        assert ok;

        var mat = new Matrix(1, 1).init(poly);
        var matI = mat.integrateIm("x");
        var matD = matI.deriveIm("x");
        System.out.println("Matrix f(x)=g`(x)" + mat);
        System.out.println("Matrix g(x)" + matI);
        System.out.println("Matrix2 f(x)" + matD);
        boolean ok2 = mat.equals(matD);
        assert ok2;
        var Lx = new Matrix(3, 3).init(
                0, 0, 0,
                0, 0, -1,
                0, 1, 0);
        var Ly = new Matrix(3, 3).init(
                0, 0, 1,
                0, 0, 0,
                -1, 0, 0);
        var Lz = new Matrix(3, 3).init(
                0, -1, 0,
                1, 0, 0,
                0, 0, 0);
        BiFunction<Matrix, Matrix, Matrix> lieIJK = (a, b) ->
                a.multiplyIm(b).minusIm(b.multiplyIm(a));
        var xx = lieIJK.apply(Ly, Lz);
        var yy = lieIJK.apply(Lz, Lx);
        var zz = lieIJK.apply(Lx, Ly);
        System.out.println("L x" + Lx);
        System.out.println("L y" + Ly);
        System.out.println("L z" + Lz);
        System.out.println("L[y,z] = yz-zy -> " + xx);
        System.out.println("L[z,x] = zx-xz -> " + yy);
        System.out.println("L[x,y] = xy-yx -> " + zz);
        boolean xOk = xx.equals(Lx);
        boolean yOk = yy.equals(Ly);
        boolean zOk = zz.equals(Lz);
        assert xOk;
        assert yOk;
        assert zOk;
    }

    static void polynomialDivide() {
        System.out.println("-- polynomialDivide --");
        var a = MVPolynomial.parse("x x + 2x + 3");
        var b = MVPolynomial.parse("x + 2");
        var ab = a.multiplyIm(b);
        var divideResult = ab.divideIm(b);
        System.out.println("a = " + a);
        System.out.println("b = " + b);
        System.out.println("ab = a * b = " + ab);
        System.out.println("ab / b = " + divideResult);
    }

    static void main(String[] args) {
        polynomialDivide();
        Matrix rot3 = calcRot3x3Matrix();
        diff();
        runCalcRot3x3Matrix();
        System.out.println(rot3);
        tstEulerMatrixPow();
        tstEulerRotMatrixTransposeProdIsIdentity();
        Matrix.setLogConsumer(System.out::println);
        testMisc();
        Matrix.printMatrixRingBufAndClear();
        //todo codegen?
        //?Ask DeepSeek.com: "Show an example of using Baker-Campbell-Hausdorff formula for composing rotation vectors"
    }

}
