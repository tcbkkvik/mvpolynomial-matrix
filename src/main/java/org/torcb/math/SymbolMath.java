package org.torcb.math;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public interface SymbolMath {
    DecimalFormat DF = new DecimalFormat("#.###", decSep());
    static DecimalFormatSymbols decSep() {
        var ds = DecimalFormatSymbols.getInstance();
        ds.setDecimalSeparator('.');
        return ds;
    }
    static boolean zero(double d) {
        return Math.abs(d) < 1e-10;
    }

    class Term implements Comparable<Term> {
        public static final Pattern EXPECT_AZ_SYMBOL_PATTERN = Pattern.compile("[a-zA-Z_ ]+.*");
        private final List<String> product = new ArrayList<>();

        public Term() {}

        public Term(String symbString) {
            if (symbString != null && !symbString.isEmpty()) {
                if (!EXPECT_AZ_SYMBOL_PATTERN.matcher(symbString).matches()) {
                    throw new IllegalArgumentException("only non-numeric symbols allowed (no scalar factor)");
                }
                build(symbString.split("[ *]"));
            }
        }

        public void build(String... arr) {
            for (var s : arr) {
                var tr = s.trim();
                if (!tr.isEmpty()) {
                    product.add(tr);
                }
            }
            sort();
        }

        public Term multiplyIm(Term o) {
            if (o == null || o.product.isEmpty()) {
                return this; // null <=> empty <=> "1"
            }
            var t = new Term();
            t.product.addAll(product);
            t.product.addAll(o.product);
            t.sort();
            return t;
        }

        public Term divideIm(Term o) {
            if (o == null) return null;
            if (o.product.isEmpty()) return this;
            var t = new Term();
            t.product.addAll(product);
            for (String s : o.product) {
                if (!t.product.remove(s)) return null;
            }
            t.sort();
            return t;
        }

        void sort() {Collections.sort(product);}

        @Override
        public final boolean equals(Object o) {
            if (!(o instanceof Term term)) return false;
            sort();
            term.sort();
            return product.equals(term.product);
        }

        @Override
        public int hashCode() {
            return product.size();
        }

        @Override
        public int compareTo(Term o) {
            int diff = product.size() - o.product.size();
            if (diff != 0) return -diff;
            sort();
            o.sort();
            var a = product.iterator();
            var b = o.product.iterator();
            while (a.hasNext() && b.hasNext()) {
                int c = a.next().compareTo(b.next());
                if (c != 0) return -c;
            }
            return 0;
        }

        @Override
        public String toString() {
            return String.join("*", product);
        }
    }

    /**
     * Multivariate Polynomial
     */
    class MVPolynomial {
        private final Map<Term, Double> map = new TreeMap<>();
        static final Pattern SCALAR_TERM_PATTERN = Pattern.compile("([^a-zA-Z_]*)(.*)");
        static final Pattern SIGN_PATTERN = Pattern.compile("[+-]");
        String label;

        public static MVPolynomial parse(String expression) {
            var sumOfTerms = new MVPolynomial();
            if (expression == null || expression.isEmpty()) return sumOfTerms;
            String[] parts = SIGN_PATTERN
                    .splitWithDelimiters(expression.replace("*", " "), 100);
            double sign = 1;
            for (String part : parts) {
                part = part.trim();
                switch (part) {
                    case "", "+" -> {}
                    case "-" -> sign *= -1;
                    default -> {
                        Matcher matcher = SCALAR_TERM_PATTERN.matcher(part);
                        if (!matcher.matches()) {
                            continue;
                        }
                        String num = matcher.group(1).trim();
                        double scalar = !num.isEmpty() ? Double.parseDouble(num) : 1;
                        sumOfTerms.add(new Term(matcher.group(2).trim()), sign * scalar);
                        sign = 1;
                    }
                }
            }
            return sumOfTerms;
        }

        public double scalarSum() {
            double sum = 0;
            for (double val : map.values()) {
                sum += val;
            }
            return sum;
        }

        public int approxSize() {
            return map.keySet().stream()
                      .map(t -> 1 + t.product.size())
                      .reduce(Integer::sum).orElse(0);
        }

        public boolean isZero() {
            removeEmpty();
            return map.isEmpty();
        }

        public MVPolynomial copy() {
            var s = new MVPolynomial();
            s.map.putAll(map);
            return s;
        }

        public MVPolynomial add(double scalar) {
            return add(new Term(), scalar);
        }

        public MVPolynomial add(String expr) {
            return add(MVPolynomial.parse(expr));
        }

        public MVPolynomial add(double scalar, String expr) {
            return add(MVPolynomial.parse(expr), scalar);
        }

        public MVPolynomial add(MVPolynomial st) {
            return add(st, 1);
        }

        public MVPolynomial minus(MVPolynomial st) {
            return add(st, -1);
        }

        public MVPolynomial add(MVPolynomial st, double scalar) {
            if (st != null)
                st.map.forEach((t, s) -> add(t, s * scalar));
            return this;
        }

        public MVPolynomial add(Term term, double scalar) {
            if (!zero(scalar) && term != null) {
                map.compute(term, (t, s) -> {
                    if (s == null) return scalar;
                    double k = s + scalar;
                    return zero(k) ? null : k;
                });
            }
            return this;
        }

        public MVPolynomial addIm(MVPolynomial other, double scalar) {
            return new MVPolynomial().add(this).add(other, scalar);
        }

        public MVPolynomial addIm(MVPolynomial other) {
            return addIm(other, 1);
        }

        public MVPolynomial minusIm(MVPolynomial other) {
            return addIm(other, -1);
        }

        public MVPolynomial negateIm() {
            return multiplyIm(-1);
        }

        public MVPolynomial multiplyIm(double scalar) {
            return multiplyIm(null, scalar);
        }

        public MVPolynomial multiplyIm(Term term, double scalar) {
            var st = new MVPolynomial();
            if (zero(scalar)) return st;
            map.forEach((t, s) ->
                    st.map.put(t.multiplyIm(term), s * scalar)
            );
            st.removeEmpty();
            return st;
        }

        public MVPolynomial multiplyIm(String expr) {
            return multiplyIm(MVPolynomial.parse(expr));
        }

        /**
         * Multiply, immutable
         */
        public MVPolynomial multiplyIm(MVPolynomial other) {
            var res = new MVPolynomial(); //empty <=> 0
            if (other == null) return res; //null <=> 0
            for (var e1 : map.entrySet()) {
                for (var e2 : other.map.entrySet()) {
                    res.add(e1.getKey().multiplyIm(e2.getKey()), e1.getValue() * e2.getValue()
                    );
                }
            }
            res.removeEmpty();
            return res;
        }

        private record STerm(double scalar, Term term) {
            STerm divide(STerm other) {
                Term res = term.divideIm(other.term);
                return res == null || zero(other.scalar()) ? null
                        : new STerm(scalar / other.scalar(), res);
            }
        }

        private static STerm highestDegreeTerm(MVPolynomial s) {
            return s == null ? null :
                    s.map.entrySet().stream().findFirst()
                         .map(e -> new STerm(e.getValue(), e.getKey()))
                         .orElse(null);
        }

        public record Divided(MVPolynomial ans, MVPolynomial remain) {
        }

        public Divided divideIm(MVPolynomial div) {
            var ans = new MVPolynomial();
            var remain = this;
            STerm d = highestDegreeTerm(div);
            if (d == null) return new Divided(ans, remain);
            while (!remain.isZero()) {
                STerm r = highestDegreeTerm(remain);
                if (r == null) break;
                STerm partial = r.divide(d);
                if (partial == null) break;
                ans.add(partial.term, partial.scalar);
                remain = remain.minusIm(div.multiplyIm(partial.term, partial.scalar));
            }
            return new Divided(ans, remain);
        }

        void removeEmpty() {
            var toRemove = map.entrySet().stream()
                              .filter(e -> zero(e.getValue()))
                              .map(Map.Entry::getKey).toList();
            toRemove.forEach(map::remove);
        }

        public MVPolynomial substituteTermsIm(String fromTerm, String toExpression) {
            return substituteTermsIm(new Term(fromTerm), MVPolynomial.parse(toExpression));
        }

        public MVPolynomial substituteTermsIm(Term sub, MVPolynomial repl) {
            var out = new MVPolynomial();
            var replCount = new AtomicInteger();
            map.forEach((t, s) -> {
                var r = t.divideIm(sub);
                if (r != null) {
                    replCount.incrementAndGet();
                    out.add(repl.multiplyIm(r, s));
                } else {
                    out.add(t, s);
                }
            });
            if (replCount.get() > 0) {
                out.label = "    //Substituted: " + sub.toString() + " -> " + repl.toString();
                out.removeEmpty();
                return out;
            }
            return this;
        }

        public MVPolynomial deriveIm(String variable) {
            var out = new MVPolynomial();
            for (Map.Entry<Term, Double> entry : map.entrySet()) {
                var t2 = new Term();
                int count = 0;
                for (var v : entry.getKey().product) {
                    if (!variable.equals(v) || ++count > 1)
                        t2.product.add(v);
                }
                out.add(t2, entry.getValue() * count);
            }
            return out;
        }

        public MVPolynomial integrateIm(String variable) {
            var out = new MVPolynomial();
            for (Map.Entry<Term, Double> entry : map.entrySet()) {
                var t2 = new Term();
                t2.product.add(variable);
                int count = 1;
                for (var v : entry.getKey().product) {
                    if (variable.equals(v)) {
                        ++count;
                    }
                    t2.product.add(v);
                }
                t2.sort();
                out.add(t2, entry.getValue() / count);
            }
            return out;
        }

        @Override
        public String toString() {
            var sb = new StringBuilder();
            map.forEach((term, value) -> {
                String v = DF.format(Math.abs(value));
                String s = term.toString();
                sb.append(sb.isEmpty()
                          ? (value < 0 ? "−" : "")
                          : (value < 0 ? " − " : " + "))
                  .append(s.isEmpty() || !"1".equals(v) ? v : "")
                  .append(s);
            });
            return sb.isEmpty() ? "0" : sb.toString();
        }

        @Override
        public final boolean equals(Object o) {
            if (o == null) {
                return isZero();
            }
            if (!(o instanceof MVPolynomial that)) {
                return false;
            }
            removeEmpty();
            that.removeEmpty();
            if (map.size() != that.map.size()) return false;
            for (var e : map.entrySet()) {
                var a = e.getValue();
                var b = that.map.get(e.getKey());
                if (a != null && b != null) {
                    if (!zero(a - b)) return false;
                } else if ((a == null) != (b == null)) {
                    return false;
                }
            }
            return true;
        }

        @Override
        public int hashCode() {
            return map.hashCode();
        }
    }

    record SubstituteTerm(Term fromTerm, MVPolynomial toExpression) {
        public static SubstituteTerm parse(String fromTerm, String toExpression) {
            return new SubstituteTerm(new Term(fromTerm), MVPolynomial.parse(toExpression));
        }

        @Override
        public String toString() {
            return "'" + fromTerm + "'->'" + toExpression + "'";
        }
    }

    class SubstituteTerms {
        public final List<SubstituteTerm> list = new ArrayList<>();

        public SubstituteTerms add(String fromTerm, String toExpression) {
            list.add(SubstituteTerm.parse(fromTerm, toExpression));
            return this;
        }
    }

    static MVPolynomial dotProduct(MVPolynomial[] a, MVPolynomial... b) {
        var res = new MVPolynomial();
        if (a.length != b.length) throw new IllegalArgumentException("different vector lengths");
        for (int i = 0; i < a.length; i++) {
            res.add(a[i].multiplyIm(b[i]));
        }
        res.removeEmpty();
        return res;
    }

    class Matrix {
        public final int nRows, nCols;
        public final MVPolynomial[] cells;
        //Equivalent semantics: cells[i]=null <=> cells[i].isZero()
        //Meta fields:
        private static final AtomicInteger ID_GENERATOR = new AtomicInteger();
        public final int id = ID_GENERATOR.incrementAndGet();
        static final LinkedList<Matrix> logRingBuf = new LinkedList<>();
        private String label = "";
        private String opArgs;
        private static volatile Consumer<String> logConsumer;

        public static void setLogConsumer(Consumer<String> c) {
            logConsumer = c;
        }

        private void logOp(String op_args) {
            opArgs = "Op: " + op_args;
            if (logConsumer != null) {
                while (logRingBuf.size() > 2) {
                    logConsumer.accept(logRingBuf.removeFirst().str());
                }
            }
        }

        public Matrix(int nRows, int nCols) {
            this.nRows = nRows;
            this.nCols = nCols;
            cells = new MVPolynomial[nRows * nCols];
            logRingBuf.add(this);
            while (logRingBuf.size() > 100) {
                logRingBuf.removeFirst();
            }
        }

        public Matrix init(MVPolynomial... arr) {
            validLength(arr.length);
            System.arraycopy(arr, 0, cells, 0, arr.length);
            return this;
        }

        public Matrix init(String... arr) {
            validLength(arr.length);
            for (int j = 0; j < arr.length; j++) {
                cells[j] = MVPolynomial.parse(arr[j]);
            }
            return this;
        }

        public Matrix init(double... values) {
            validLength(values.length);
            for (int j = 0; j < values.length; j++) {
                cells[j] = new MVPolynomial().add(values[j]);
            }
            return this;
        }

        private void validLength(int length) {
            if (length != nRows * nCols)
                throw new IllegalArgumentException(
                        String.format("Expected nRows*nCols %d*%d=%d, Got:%d"
                                , nRows, nCols, cells.length, length));
        }

        public static Matrix identity(int dim) {
            return diagonal(dim, new MVPolynomial().add(1));
        }

        public static Matrix diagonal(int dim, MVPolynomial value) {
            var out = new Matrix(dim, dim);
            for (int pos = 0; pos < out.cells.length; pos += dim + 1) {
                out.cells[pos] = value;
            }
            out.logOp("diagonal(dim " + dim + ", val '" + value + "')");
            return out;
        }

        public static Matrix init3x3(
                String a, String b, String c,
                String d, String e, String f,
                String g, String h, String i) {
            return new Matrix(3, 3)
                    .init(a, b, c, d, e, f, g, h, i);
        }

        public Matrix label(String label) {
            if (label != null) this.label = label;
            return this;
        }

        public Matrix transposeIm() {
            var out = new Matrix(nCols, nRows);
            iterateNonNull((pos, row, col, cell) ->
                    out.cells[col * nRows + row] = cell.copy());
            if (!isEmpty(label))
                out.label("tr(" + label + ")");
            out.logOp(id + ".transpose()");
            return out;
        }

        public Matrix minusIm(Matrix right) {
            return addIm(right, -1);
        }

        public Matrix addIm(Matrix right) {
            return addIm(right, 1);
        }

        public Matrix addIm(Matrix right, double scalar) {
            if (nRows != right.nRows || nCols != right.nCols)
                throw new IllegalArgumentException("Dimension mismatch");
            var out = new Matrix(nRows, nCols);
            iterate((pos, r, c, a) -> {
                var b = right.cells[pos];
                out.cells[pos] = (a == null && b == null) ?
                        null : new MVPolynomial().add(a).add(b, scalar);
            });
            out.logOp(id + ".add(" + scalar + " * matrix " + right.id + ")");
            return out;
        }

        public Matrix multiplyIm(double scalar) {
            return multiplyIm(scalar, null);
        }

        /**
         * @param scalar multiplier
         * @param term   if null => scalar op only
         * @return new matrix
         */
        public Matrix multiplyIm(double scalar, Term term) {
            var out = new Matrix(nRows, nCols);
            iterateNonNull((pos, row, col, cell) ->
                    out.cells[pos] = cell.multiplyIm(term, scalar));
            var t = term == null ? "" : " * " + term;
            out.logOp(id + ".multiplyIm(" + scalar + t + ")");
            return out;
        }

        public Matrix multiplyIm(String expression) {
            return multiplyIm(MVPolynomial.parse(expression));
        }

        public Matrix multiplyIm(MVPolynomial expression) {
            var out = new Matrix(nRows, nCols);
            if (expression == null || expression.isZero()) return out;
            iterateNonNull((pos, row, col, cell) ->
                    out.cells[pos] = cell.multiplyIm(expression));
            out.logOp(id + ".multiply(expr '" + expression + "')");
            return out;
        }

        public Matrix multiplyIm(Matrix right) {
            if (nCols != right.nRows) {
                throw new IllegalArgumentException("multiply: nCols != other.nRows");
            }
            var out = new Matrix(nRows, right.nCols);
            out.iterate((pos, row, col, cell) -> {
                var elem = out.cells[pos] = new MVPolynomial();
                int rightPos = col;
                int leftPos = row * nCols;
                for (int i = 0; i < nCols; i++) {
                    var l = cells[leftPos++];
                    var r = right.cells[rightPos];
                    if (l != null && r != null) { //null <=> 0
                        elem.add(l.multiplyIm(r));
                    }
                    rightPos += out.nCols;
                }
            });
            out.logOp(id + ".multiply(matrix " + right.id + ")");
            return out;
        }

        public Matrix substituteTermsIm(SubstituteTerms subst) {
            Matrix mat = this;
            for (var st : subst.list) {
                mat = mat.substituteTermsIm(st);
            }
            mat.label(label);
            return mat;
        }

        public Matrix substituteTermsIm(String fromTerm, String toExpression) {
            return substituteTermsIm(SubstituteTerm.parse(fromTerm, toExpression));
        }

        public Matrix substituteTermsIm(SubstituteTerm st) {
            var out = new Matrix(nRows, nCols);
            var replCount = new AtomicInteger();
            iterateNonNull((pos, row, col, cell) -> {
                var replaced = cell.substituteTermsIm(st.fromTerm, st.toExpression);
                if (replaced.approxSize() <= cell.approxSize()) {
                    out.cells[pos] = replaced;
                    replCount.incrementAndGet();
                } else {
                    out.cells[pos] = cell;
                }
            });
            if (replCount.get() > 0) {
                out.label(label);
                out.logOp(id + ".replace(" + st + ") #:" + replCount.get());
                return out;
            }
            ID_GENERATOR.decrementAndGet();
            logRingBuf.removeLast();
            return this;
        }

        public Matrix deriveIm(String variable) {
            var out = new Matrix(nRows, nCols);
            iterateNonNull((pos, row, col, cell) ->
                    out.cells[pos] = cell.deriveIm(variable));
            out.logOp(id + ".derive(" + variable + ")");
            return out;
        }

        public Matrix integrateIm(String variable) {
            var out = new Matrix(nRows, nCols);
            iterateNonNull((pos, row, col, cell) ->
                    out.cells[pos] = cell.integrateIm(variable));
            out.logOp(id + ".integrate(" + variable + ")");
            return out;
        }

        public double[] getAllScalars() {
            double[] scalars = new double[cells.length];
            for (int l = 0; l < cells.length; l++) {
                scalars[l] = cells[l].scalarSum();
            }
            return scalars;
        }

        public <T> T op(Function<Matrix, T> op) {
            return op.apply(this);
        }

        @Override
        public final boolean equals(Object o) {
            if (!(o instanceof Matrix matrix)) return false;
            if (nRows != matrix.nRows || nCols != matrix.nCols) {
                return false;
            }
            for (int pos = 0; pos < cells.length; pos++) {
                MVPolynomial a = cells[pos];
                MVPolynomial b = matrix.cells[pos];
                if (a != null) {
                    if (!a.equals(b)) return false;
                } else {
                    if (b != null && !b.isZero()) return false;
                }
            }
            return true;
        }

        @Override
        public int hashCode() {
            return id;
        }

        public String str() {
            var s = isEmpty(opArgs) ? "" : "\n " + opArgs + " ->\n";
            return s + this;
        }

        @Override
        public String toString() {
            var lb = isEmpty(label) ? "" : " " + label;
            var sb = new StringBuilder(" Matrix{" + id + lb + "}\n  ");
            iterate((pos, row, col, cell) -> {
                if (pos > 0 && col == 0) sb.append("\n  ");
                sb.append(cell == null ? "0" : cell)
                  .append(col < nCols - 1 ? ",  " : ";");
            });
            return sb.toString();
        }

        public static void printMatrixRingBufAndClear() {
            while (!logRingBuf.isEmpty()) {
                String s = logRingBuf.removeFirst().str();
                if (logConsumer != null) logConsumer.accept(s);
                else System.out.println(s);
            }
        }

        public void iterateNonNull(RowCol action) {
            iterate(action, false);
        }

        public void iterate(RowCol action) {
            iterate(action, true);
        }

        public void iterate(RowCol action, boolean allowNull) {
            int row = 0, col = 0;
            for (int pos = 0; pos < cells.length; pos++) {
                if (cells[pos] != null || allowNull) {
                    action.accept(pos, row, col, cells[pos]);
                }
                if (++col < nCols) continue;
                ++row;
                col = 0;
            }
        }

        public MVPolynomial cell(int row, int col) {
            if (col < 0 || col >= nCols)
                throw new IllegalArgumentException("Col outside [0.." + nCols + ")");
            return cells[row * nCols + col];
        }

        @FunctionalInterface
        public interface RowCol {
            void accept(int pos, int row, int col, MVPolynomial cell);
        }

        /**
         * Get diagonal sum
         *
         * @param matrix source (immutable)
         * @return SumOfTerms
         */
        public static MVPolynomial trace(Matrix matrix) {
            var st = new MVPolynomial();
            if (matrix == null) return st;
            int N = Math.min(matrix.nRows, matrix.nCols);
            for (int i = 0; i < N; i++) {
                st.add(matrix.cell(i, i));
            }
            //From rotate3dMatrix R: angle fi = arcCos( (trace(R)-1)/2 )
            return st;
        }
    }

    static boolean isEmpty(String s) {
        return s == null || s.isEmpty();
    }

}
