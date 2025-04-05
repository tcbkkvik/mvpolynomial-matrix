
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
        var replace_ii = new SubstituteTerms()
                .add("i i", "1 - j j - k k")
                .add("j j", "1 - i i - k k");
        var replace_cos2 = new SubstituteTerms()
                .add("cos cos", "1 - sin sin");
        var L = Matrix.init3x3(
                "0", "-k", "j",
                "k", "0", "-i",
                "-j", "i", "0").label("L (CrossProduct matrix; Lie algebra)");
        var L_sin = L.multiplyIm("sin").label("L*sin");
        var LL = L.multiplyIm(L).label("L*L").substituteTermsIm(replace_ii);
        var LL_1_cos = LL.multiplyIm("1 - cos").label("LL*(1-cos)");// 1 - cos;
        Matrix rotateM = Matrix.identity(3).addIm(L_sin).addIm(LL_1_cos)
                               .label("R = I + L*sin - L*L*(1-cos)  //Rodrigues formula");
        // - R*tr(T)==I:
        var idRot = rotateM.multiplyIm(rotateM.transposeIm()).label("I = R*tr(R)");
        var idRot2 = idRot.substituteTermsIm(replace_ii);
        var idRot3 = idRot2.substituteTermsIm(replace_cos2);
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
 Matrix{3 L*L}
  −k*k − j*j,  i*j,  i*k;
  i*j,  −k*k − i*i,  j*k;
  i*k,  j*k,  −j*j − i*i;

 Op: 3.replace('i*i'->'−k*k − j*j + 1') #:9 ->
 Matrix{4 L*L}
  −k*k − j*j,  i*j,  i*k;
  i*j,  j*j − 1,  j*k;
  i*k,  j*k,  k*k − 1;

 Op: 4.replace('j*j'->'−k*k − i*i + 1') #:8 ->
 Matrix{5 L*L}
  i*i − 1,  i*j,  i*k;
  i*j,  j*j − 1,  j*k;
  i*k,  j*k,  k*k − 1;

 Op: 5.multiply(expr '−cos + 1') ->
 Matrix{6 LL*(1-cos)}
  −cos*i*i + i*i + cos − 1,  −cos*i*j + i*j,  −cos*i*k + i*k;
  −cos*i*j + i*j,  −cos*j*j + j*j + cos − 1,  −cos*j*k + j*k;
  −cos*i*k + i*k,  −cos*j*k + j*k,  −cos*k*k + k*k + cos − 1;

 Op: diagonal(dim 3, val '1') ->
 Matrix{7}
  1,  0,  0;
  0,  1,  0;
  0,  0,  1;

 Op: 7.add(1.0 * matrix 2) ->
 Matrix{8}
  1,  −k*sin,  j*sin;
  k*sin,  1,  −i*sin;
  −j*sin,  i*sin,  1;

 Op: 8.add(1.0 * matrix 6) ->
 Matrix{9 Rotate3D = I + L*sin - L*L*(1-cos)  //Rodrigues formula}
  −cos*i*i + i*i + cos,  −cos*i*j − k*sin + i*j,  −cos*i*k + j*sin + i*k;
  −cos*i*j + k*sin + i*j,  −cos*j*j + j*j + cos,  −cos*j*k + j*k − i*sin;
  −cos*i*k − j*sin + i*k,  −cos*j*k + j*k + i*sin,  −cos*k*k + k*k + cos;

 Op: 9.transpose() ->
 Matrix{10 tr(Rotate3D = I + L*sin - L*L*(1-cos)  //Rodrigues formula)}
  −cos*i*i + i*i + cos,  −cos*i*j + k*sin + i*j,  −cos*i*k − j*sin + i*k;
  −cos*i*j − k*sin + i*j,  −cos*j*j + j*j + cos,  −cos*j*k + j*k + i*sin;
  −cos*i*k + j*sin + i*k,  −cos*j*k + j*k − i*sin,  −cos*k*k + k*k + cos;

 Op: 9.multiply(matrix 10) ->
 Matrix{11 I = Rotate3D * transpose(Rotate3D)}
  cos*cos*i*i*k*k + cos*cos*i*i*j*j + cos*cos*i*i*i*i − 2cos*i*i*k*k − 2cos*i*i*j*j − 2cos*i*i*i*i + k*k*sin*sin + j*j*sin*sin + i*i*k*k + i*i*j*j + i*i*i*i − 2cos*cos*i*i + 2cos*i*i + cos*cos,  cos*cos*i*j*k*k + cos*cos*i*j*j*j + cos*cos*i*i*i*j − 2cos*i*j*k*k − 2cos*i*j*j*j − 2cos*i*i*i*j − i*j*sin*sin + i*j*k*k + i*j*j*j + i*i*i*j − 2cos*cos*i*j + 2cos*i*j,  cos*cos*i*k*k*k + cos*cos*i*j*j*k + cos*cos*i*i*i*k − 2cos*i*k*k*k − 2cos*i*j*j*k − 2cos*i*i*i*k − i*k*sin*sin + i*k*k*k + i*j*j*k + i*i*i*k − 2cos*cos*i*k + 2cos*i*k;
  cos*cos*i*j*k*k + cos*cos*i*j*j*j + cos*cos*i*i*i*j − 2cos*i*j*k*k − 2cos*i*j*j*j − 2cos*i*i*i*j − i*j*sin*sin + i*j*k*k + i*j*j*j + i*i*i*j − 2cos*cos*i*j + 2cos*i*j,  cos*cos*j*j*k*k + cos*cos*j*j*j*j + cos*cos*i*i*j*j − 2cos*j*j*k*k − 2cos*j*j*j*j − 2cos*i*i*j*j + k*k*sin*sin + j*j*k*k + j*j*j*j + i*i*sin*sin + i*i*j*j − 2cos*cos*j*j + 2cos*j*j + cos*cos,  cos*cos*j*k*k*k + cos*cos*j*j*j*k + cos*cos*i*i*j*k − 2cos*j*k*k*k − 2cos*j*j*j*k − 2cos*i*i*j*k − j*k*sin*sin + j*k*k*k + j*j*j*k + i*i*j*k − 2cos*cos*j*k + 2cos*j*k;
  cos*cos*i*k*k*k + cos*cos*i*j*j*k + cos*cos*i*i*i*k − 2cos*i*k*k*k − 2cos*i*j*j*k − 2cos*i*i*i*k − i*k*sin*sin + i*k*k*k + i*j*j*k + i*i*i*k − 2cos*cos*i*k + 2cos*i*k,  cos*cos*j*k*k*k + cos*cos*j*j*j*k + cos*cos*i*i*j*k − 2cos*j*k*k*k − 2cos*j*j*j*k − 2cos*i*i*j*k − j*k*sin*sin + j*k*k*k + j*j*j*k + i*i*j*k − 2cos*cos*j*k + 2cos*j*k,  cos*cos*k*k*k*k + cos*cos*j*j*k*k + cos*cos*i*i*k*k − 2cos*k*k*k*k − 2cos*j*j*k*k − 2cos*i*i*k*k + k*k*k*k + j*j*sin*sin + j*j*k*k + i*i*sin*sin + i*i*k*k − 2cos*cos*k*k + 2cos*k*k + cos*cos;

 Op: 11.replace('i*i'->'−k*k − j*j + 1') #:8 ->
 Matrix{12 I = Rotate3D * transpose(Rotate3D)}
  cos*cos*i*i*k*k + cos*cos*i*i*j*j + cos*cos*i*i*i*i − 2cos*i*i*k*k − 2cos*i*i*j*j − 2cos*i*i*i*i + k*k*sin*sin + j*j*sin*sin + i*i*k*k + i*i*j*j + i*i*i*i − 2cos*cos*i*i + 2cos*i*i + cos*cos,  −i*j*sin*sin − cos*cos*i*j + i*j,  −i*k*sin*sin − cos*cos*i*k + i*k;
  −i*j*sin*sin − cos*cos*i*j + i*j,  −j*j*sin*sin − cos*cos*j*j + sin*sin + j*j + cos*cos,  −j*k*sin*sin − cos*cos*j*k + j*k;
  −i*k*sin*sin − cos*cos*i*k + i*k,  −j*k*sin*sin − cos*cos*j*k + j*k,  −k*k*sin*sin − cos*cos*k*k + sin*sin + k*k + cos*cos;

 Op: 12.replace('j*j'->'−k*k − i*i + 1') #:8 ->
 Matrix{13 I = Rotate3D * transpose(Rotate3D)}
  −i*i*sin*sin − cos*cos*i*i + sin*sin + i*i + cos*cos,  −i*j*sin*sin − cos*cos*i*j + i*j,  −i*k*sin*sin − cos*cos*i*k + i*k;
  −i*j*sin*sin − cos*cos*i*j + i*j,  −j*j*sin*sin − cos*cos*j*j + sin*sin + j*j + cos*cos,  −j*k*sin*sin − cos*cos*j*k + j*k;
  −i*k*sin*sin − cos*cos*i*k + i*k,  −j*k*sin*sin − cos*cos*j*k + j*k,  −k*k*sin*sin − cos*cos*k*k + sin*sin + k*k + cos*cos;

 Op: 13.replace('cos*cos'->'−sin*sin + 1') #:9 ->
 Matrix{14 I = Rotate3D * transpose(Rotate3D)}
  1,  0,  0;
  0,  1,  0;
  0,  0,  1;
````

Initiated via:
mvn archetype:generate -DgroupId=org.torcb.project -DartifactId=__artifact__ -DarchetypeArtifactId=maven-archetype-quickstart -DinteractiveMode=false
