package org.torcb.math;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
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
        SUBSTITUTE_RULES.set(new SubstituteTerms()
                .add("j j", "1 - i i - k k")
                .add("cos cos", "1 - sin sin"));

        var L = Matrix.init3x3(
                "0", "-k", "j",
                "k", "0", "-i",
                "-j", "i", "0").label("L (CrossProd matrix)");
        var L_sin = L.multiplyIm("sin").label("L*sin");
        var LL = L.multiplyIm(L).label("L*L");
        var LL_1_cos = LL.multiplyIm("1 - cos").label("LL*(1-cos)");// 1 - cos;
        Matrix identity3D = Matrix.identity(3);
        Matrix rotateM = identity3D.addIm(L_sin).addIm(LL_1_cos)
                                   .label("Rotate = I + L*sin - L*L*(1-cos)  //Rodrigues formula");
        Matrix idRot = rotateM.multiplyIm(rotateM.transposeIm()).label("I = Rotate * transpose(Rotate)");
        MVPolynomial det = rotateM.determinant();

        assertEquals(identity3D, idRot);
        assertEquals(new MVPolynomial().add(1), det);
    }

}
