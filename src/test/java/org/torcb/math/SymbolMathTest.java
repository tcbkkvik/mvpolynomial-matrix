package org.torcb.math;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;
import static org.torcb.math.SymbolMath.*;

public class SymbolMathTest {

    @Test
    void testTerm() {
        var b = new Term("b");
        var a = new Term("a");
        Term ab = a.multiplyIm(b);
        var expect_a = ab.divideIm(b);
        assertEquals(new Term("a * b"), ab);
        assertEquals(expect_a, a);
    }

    @Test
    void testSumOfTerms() {
        var s1 = MVPolynomial.parse("a b + 2 c");
        var s2 = MVPolynomial.parse("1 - cos cos")
                             .substituteTermsIm("cos cos", "1 - sin sin");
        var s12 = s1.multiplyIm(s2);
        var s12exp = MVPolynomial.parse("a*b*sin*sin + 2c*sin*sin");
        assertEquals(s12exp, s12);
        assertTrue(s12exp.substituteTermsIm("a b", "-2c").isZero());

        var a = MVPolynomial.parse("x x + 2x + 3");
        var b = MVPolynomial.parse("x + 2");
        assertEquals(a, a.minus(b).add(b));
        var a2 = a.addIm(a);
        assertEquals(a.multiplyIm(2), a2);
        var ab = a.multiplyIm(b);
        var ab_b = ab.divideIm(b).ans();
        assertEquals(a, ab_b);

        var poly = MVPolynomial.parse("3 x x a + 2 x b+ 9 c");
        var expect = MVPolynomial.parse("2b + 6a*x");
        var derived = poly.deriveIm("x");
        assertEquals(expect, derived);

        var mat = new Matrix(1, 1).init(poly);
        var matI = mat.integrateIm("x");
        var matD = matI.deriveIm("x");
        assertEquals(mat, matD);
    }

    @Test
    void testSumOfTerms2() {
        /* https://thenumb.at/Exponential-Rotations/#euler-angles
                      rotMat = L*L − L*L*cos + I + L*sin
            transpose(rotMat)= L*L − L*L*cos + I − L*sin
            transpose(rotMat) * rotMat
             = I*I + L*L*L*L + 2I*L*L − L*L*sin*sin + L*L*L*L*cos*cos − 2L*L*L*L*cos − 2I*L*L*cos
             = −L*L + I*I − L*L*cos*cos + 2L*L*cos + 2I*L*L − L*L*sin*sin − 2I*L*L*cos
             = I*I − 2L*L + 2I*L*L + 2L*L*cos − 2I*L*L*cos
             = I*I
             = I
        */
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
        assertEquals(MVPolynomial.parse("I"), prod_5);
    }

    @Test
    void testMatrix() {
        Matrix.logRingBuf.clear();
        SubstituteRules.set(new SubstituteTerms()
                .add("i i", "1 - j j - k k")
                .add("j j", "1 - i i - k k")
                .add("cos cos", "1 - sin sin"));

        Matrix identity3d = Matrix.identity(3);
        var L = Matrix.init3x3(
                "0", "-k", "j",
                "k", "0", "-i",
                "-j", "i", "0").label("L (CrossProd matrix)");
        var L_sin = L.multiplyIm("sin").label("L*sin");
        var LL = L.multiplyIm(L).label("L*L");
        var LL_1_cos = LL.multiplyIm("1 - cos").label("LL*(1-cos)");// 1 - cos;
        Matrix rotateM = identity3d.addIm(L_sin).addIm(LL_1_cos)
                                   .label("Rotate = I + L*sin + L*L*(1-cos)  //Euler-Rodrigues formula");
        Matrix idRot = rotateM.multiplyIm(rotateM.transposeIm()).label("I = Rotate * transpose(Rotate)");

        var mat2 = Matrix.parse(rotateM.toString()).label("parsed from " + rotateM.id);
        MVPolynomial det = rotateM.determinant();
        Matrix.printMatrixRingBufAndClear();

        System.out.printf("\n Matrix %s determinant: %s\n", rotateM.id, det);
        assertEquals(rotateM, mat2);
        assertEquals(identity3d, idRot);
        assertEquals(new MVPolynomial().add(1), det);
    }

    @Test
    void testParseFull() throws NumberFormatException {
        try {
            MVPolynomial.parse("sld (( [&& sdf fj032 j-as");
        } catch (Exception ex) {
            //noinspection ThrowablePrintedToSystemOut
            System.out.println(ex);
        }
        SubstituteRules.set(new SubstituteTerms()
                .add("cos cos", "1 - sin sin"));

        MVPolynomial res = MVPolynomial.parse("6a - -4 * -a");
        assertEquals(MVPolynomial.parse("2a"), res);
        MVPolynomial poly = MVPolynomial
                .parse("1 + (a + cos)(a - cos) + (a+1)(a-1) - sin sin")
                .substituteTermsIm();

        System.out.println(poly);
        var expect = MVPolynomial.parse("2 a*a  -1");
        assertEquals(expect, poly);
    }
}
