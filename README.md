
# Multivariate Polynomial Matrix - Basic Algebra on Symbolic Expressions and Matrices


## Features
* Basic Algebraic operations on Symbolic Expressions(+-*/) and Matrices(+-*transpose) 
* Expressions normalized as sum of weighted Terms:  
  k1\*t1 + k2\*t2 + ... + ki\*ti + ... + kN\*tN  
  ,  double ki, Term ti,  class Term{List\<String> symbols}
* Simplify by Term substitution (Replace Terms by Expressions: simple, deterministic, generic, powerful)
 

## Example - Simplify expression
````
     // https://thenumb.at/Exponential-Rotations/#euler-angles
        var rotMat = SumOfTerms
                .parse("1 - cos")
                .multiplyIm("L L")
                .add("I + L sin");
        var rotMatTranspose = SumOfTerms
                .parse("1 - cos")
                .multiplyIm(" L L")
                .add("I - L sin");
        var prod_1 = rotMatTranspose.multiplyIm(rotMat);
        var prod_2 = prod_1.substituteTermsIm("L L L L", "- L L");
        var prod_3 = prod_2.substituteTermsIm("sin sin", "1 - cos cos");
        var prod_4 = prod_3.substituteTermsIm("I L", "L");
        var prod_5 = prod_4.substituteTermsIm("I I", "I");
        
        //Result = Identity matrix 'I' expected
        System.out.println("tr(rotMat)= " + rotMatTranspose);
        System.out.println("rotMat = " + rotMat);
        System.out.println("tr(rotMat) * rotMat");
        System.out.println(" = " + prod_1);
        System.out.println(" = " + prod_2);
        System.out.println(" = " + prod_3);
        System.out.println(" = " + prod_4);
        System.out.println(" = " + prod_5);
````
Output:
````
tr(rotMat)= I − L*sin + L*L − L*L*cos
rotMat = I + L*sin + L*L − L*L*cos
tr(rotMat) * rotMat
 = I*I + 2I*L*L + L*L*L*L − L*L*sin*sin − 2I*L*L*cos − 2L*L*L*L*cos + L*L*L*L*cos*cos
 = −L*L + I*I + 2L*L*cos + 2I*L*L − L*L*cos*cos − 2I*L*L*cos − L*L*sin*sin
 = I*I − 2L*L + 2I*L*L + 2L*L*cos − 2I*L*L*cos
 = I*I
 = I
````

## Example - 3D rotation: I = M*transpose(M)
````
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
````
Output:
````
 Matrix{1 L (CrossProduct matrix; Lie algebra)}
  0,  −k,  j;
  k,  0,  −i;
  −j,  i,  0;

 Op: 1.multiply(expr 'sin') ->
 Matrix{2 L*sin}
  0,  −k*sin,  j*sin;
  k*sin,  0,  −i*sin;
  −j*sin,  i*sin,  0;

 Op: 1.multiply(matrix 1) ->
 Matrix{3}
  −k*k − j*j,  i*j,  i*k;
  i*j,  −k*k − i*i,  j*k;
  i*k,  j*k,  −j*j − i*i;

 Op: 3.replace('i*i'->'−k*k − j*j + 1') #:2 ->
 Matrix{4 L*L}
  −k*k − j*j,  i*j,  i*k;
  i*j,  j*j − 1,  j*k;
  i*k,  j*k,  k*k − 1;

 Op: 4.multiply(expr '−cos + 1') ->
 Matrix{5 LL*(1-cos)}
  cos*k*k + cos*j*j − k*k − j*j,  −cos*i*j + i*j,  −cos*i*k + i*k;
  −cos*i*j + i*j,  −cos*j*j + j*j + cos − 1,  −cos*j*k + j*k;
  −cos*i*k + i*k,  −cos*j*k + j*k,  −cos*k*k + k*k + cos − 1;

 Op: diagonal(dim 3, val '1') ->
 Matrix{6}
  1,  0,  0;
  0,  1,  0;
  0,  0,  1;

 Op: 6.add(1.0 * matrix 2) ->
 Matrix{7}
  1,  −k*sin,  j*sin;
  k*sin,  1,  −i*sin;
  −j*sin,  i*sin,  1;

 Op: 7.add(1.0 * matrix 5) ->
 Matrix{8 Rotate3D = I + L*sin - L*L*(1-cos)  //Rodrigues formula}
  cos*k*k + cos*j*j − k*k − j*j + 1,  −cos*i*j − k*sin + i*j,  −cos*i*k + j*sin + i*k;
  −cos*i*j + k*sin + i*j,  −cos*j*j + j*j + cos,  −cos*j*k + j*k − i*sin;
  −cos*i*k − j*sin + i*k,  −cos*j*k + j*k + i*sin,  −cos*k*k + k*k + cos;

 Op: 8.transpose() ->
 Matrix{9 tr(Rotate3D = I + L*sin - L*L*(1-cos)  //Rodrigues formula)}
  cos*k*k + cos*j*j − k*k − j*j + 1,  −cos*i*j + k*sin + i*j,  −cos*i*k − j*sin + i*k;
  −cos*i*j − k*sin + i*j,  −cos*j*j + j*j + cos,  −cos*j*k + j*k + i*sin;
  −cos*i*k + j*sin + i*k,  −cos*j*k + j*k − i*sin,  −cos*k*k + k*k + cos;

 Op: 8.multiply(matrix 9) ->
 Matrix{10}
  cos*cos*k*k*k*k + 2cos*cos*j*j*k*k + cos*cos*j*j*j*j + cos*cos*i*i*k*k + cos*cos*i*i*j*j − 2cos*k*k*k*k − 4cos*j*j*k*k − 2cos*j*j*j*j − 2cos*i*i*k*k − 2cos*i*i*j*j + k*k*sin*sin + k*k*k*k + j*j*sin*sin + 2j*j*k*k + j*j*j*j + i*i*k*k + i*i*j*j + 2cos*k*k + 2cos*j*j − 2k*k − 2j*j + 1,  cos*k*k*k*sin + cos*j*j*k*sin + cos*i*i*k*sin − k*k*k*sin − j*j*k*sin − i*j*sin*sin − i*i*k*sin − cos*cos*i*j − cos*k*sin + k*sin + i*j,  −cos*j*k*k*sin − cos*j*j*j*sin − cos*i*i*j*sin + j*k*k*sin + j*j*j*sin − i*k*sin*sin + i*i*j*sin − cos*cos*i*k + cos*j*sin − j*sin + i*k;
  cos*k*k*k*sin + cos*j*j*k*sin + cos*i*i*k*sin − k*k*k*sin − j*j*k*sin − i*j*sin*sin − i*i*k*sin − cos*cos*i*j − cos*k*sin + k*sin + i*j,  cos*cos*j*j*k*k + cos*cos*j*j*j*j + cos*cos*i*i*j*j − 2cos*j*j*k*k − 2cos*j*j*j*j − 2cos*i*i*j*j + k*k*sin*sin + j*j*k*k + j*j*j*j + i*i*sin*sin + i*i*j*j − 2cos*cos*j*j + 2cos*j*j + cos*cos,  cos*cos*j*k*k*k + cos*cos*j*j*j*k + cos*cos*i*i*j*k − 2cos*j*k*k*k − 2cos*j*j*j*k − 2cos*i*i*j*k − j*k*sin*sin + j*k*k*k + j*j*j*k + i*i*j*k − 2cos*cos*j*k + 2cos*j*k;
  −cos*j*k*k*sin − cos*j*j*j*sin − cos*i*i*j*sin + j*k*k*sin + j*j*j*sin − i*k*sin*sin + i*i*j*sin − cos*cos*i*k + cos*j*sin − j*sin + i*k,  cos*cos*j*k*k*k + cos*cos*j*j*j*k + cos*cos*i*i*j*k − 2cos*j*k*k*k − 2cos*j*j*j*k − 2cos*i*i*j*k − j*k*sin*sin + j*k*k*k + j*j*j*k + i*i*j*k − 2cos*cos*j*k + 2cos*j*k,  cos*cos*k*k*k*k + cos*cos*j*j*k*k + cos*cos*i*i*k*k − 2cos*k*k*k*k − 2cos*j*j*k*k − 2cos*i*i*k*k + k*k*k*k + j*j*sin*sin + j*j*k*k + i*i*sin*sin + i*i*k*k − 2cos*cos*k*k + 2cos*k*k + cos*cos;

 Op: 10.replace('i*i'->'−k*k − j*j + 1') #:9 ->
 Matrix{11}
  k*k*sin*sin + j*j*sin*sin + cos*cos*k*k + cos*cos*j*j − k*k − j*j + 1,  −i*j*sin*sin − cos*cos*i*j + i*j,  −i*k*sin*sin − cos*cos*i*k + i*k;
  −i*j*sin*sin − cos*cos*i*j + i*j,  −j*j*sin*sin − cos*cos*j*j + sin*sin + j*j + cos*cos,  −j*k*sin*sin − cos*cos*j*k + j*k;
  −i*k*sin*sin − cos*cos*i*k + i*k,  −j*k*sin*sin − cos*cos*j*k + j*k,  −k*k*sin*sin − cos*cos*k*k + sin*sin + k*k + cos*cos;

 Op: 11.replace('cos*cos'->'−sin*sin + 1') #:9 ->
 Matrix{12 I = Rotate3D * transpose(Rotate3D)}
  1,  0,  0;
  0,  1,  0;
  0,  0,  1;
````

Initiated via:
mvn archetype:generate -DgroupId=org.torcb.project -DartifactId=__artifact__ -DarchetypeArtifactId=maven-archetype-quickstart -DinteractiveMode=false
