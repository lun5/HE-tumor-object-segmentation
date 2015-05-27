package matrix;

import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.util.StringTokenizer;

import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.BufferedReader;

import util.Utility;



public class Matrix implements Cloneable, java.io.Serializable {

	/* ------------------------
   Class variables
	 * ------------------------ */

	/** Array for internal storage of elements.
   @serial internal array storage.
	 */
	private double[][] A;

	/** Row and column dimensions.
   @serial row dimension.
   @serial column dimension.
	 */
	private int m, n;

	int unsegmentedRegionNo;
	/* ------------------------
   Constructors
	 * ------------------------ */

	/** Construct an m-by-n matrix of zeros. 
   @param m    Number of rows.
   @param n    Number of colums.
	 */

	public Matrix (int m, int n) {
		this.m = m;
		this.n = n;
		A = new double[m][n];
	}

	/** Construct an m-by-n constant matrix.
   @param m    Number of rows.
   @param n    Number of colums.
   @param s    Fill the matrix with this scalar value.
	 */

	public Matrix (int m, int n, double s) {
		this.m = m;
		this.n = n;
		A = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = s;
			}
		}
	}

	/** Construct a matrix from a 2-D array.
   @param A    Two-dimensional array of doubles.
   @exception  IllegalArgumentException All rows must have the same length
   @see        #constructWithCopy
	 */

	public Matrix (double[][] A) {
		m = A.length;
		n = A[0].length;
		for (int i = 0; i < m; i++) {
			if (A[i].length != n) {
				throw new IllegalArgumentException("All rows must have the same length.");
			}
		}
		this.A = A;
	}

	/** Construct a matrix quickly without checking arguments.
   @param A    Two-dimensional array of doubles.
   @param m    Number of rows.
   @param n    Number of colums.
	 */

	public Matrix (double[][] A, int m, int n) {
		this.A = A;
		this.m = m;
		this.n = n;
	}

	/** Construct a matrix from a one-dimensional packed array
   @param vals One-dimensional array of doubles, packed by columns (ala Fortran).
   @param m    Number of rows.
   @exception  IllegalArgumentException Array length must be a multiple of m.
	 */

	public Matrix (double vals[], int m) {
		this.m = m;
		n = (m != 0 ? vals.length/m : 0);
		if (m*n != vals.length) {
			throw new IllegalArgumentException("Array length must be a multiple of m.");
		}
		A = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = vals[i+j*m];
			}
		}
	}

	/* ------------------------
   Public Methods
	 * ------------------------ */

	/** Construct a matrix from a copy of a 2-D array.
   @param A    Two-dimensional array of doubles.
   @exception  IllegalArgumentException All rows must have the same length
	 */

	public static Matrix constructWithCopy(double[][] A) {
		int m = A.length;
		int n = A[0].length;
		Matrix X = new Matrix(m,n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			if (A[i].length != n) {
				throw new IllegalArgumentException
				("All rows must have the same length.");
			}
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j];
			}
		}
		return X;
	}

	/** Make a deep copy of a matrix
	 */

	public Matrix copy () {
		Matrix X = new Matrix(m,n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j];
			}
		}
		X.setUnsegmentedRegionNo(this.unsegmentedRegionNo);
		return X;
	}

	/** Clone the Matrix object.
	 */

	public Object clone () {
		return this.copy();
	}

	/** Access the internal two-dimensional array.
   @return     Pointer to the two-dimensional array of matrix elements.
	 */

	public double[][] getArray () {
		return A;
	}

	/** Copy the internal two-dimensional array.
   @return     Two-dimensional array copy of matrix elements.
	 */

	public double[][] getArrayCopy () {
		double[][] C = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j];
			}
		}
		return C;
	}

	/** Make a one-dimensional column packed copy of the internal array.
   @return     Matrix elements packed in a one-dimensional array by columns.
	 */

	public double[] getColumnPackedCopy () {
		double[] vals = new double[m*n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				vals[i+j*m] = A[i][j];
			}
		}
		return vals;
	}



	public double[] getColumn(int index) {
		double[] vals = new double[n];
		for (int i = 0; i < m; i++) {
			if(i==index){
				for (int j = 0; j < n; j++) {
					vals[j] = A[i][j];
				}
				return vals;
			}
		}
		return null;
	}

	/** Make a one-dimensional row packed copy of the internal array.
   @return     Matrix elements packed in a one-dimensional array by rows.
	 */

	public double[] getRowPackedCopy () {
		double[] vals = new double[m*n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				vals[i*n+j] = A[i][j];
			}
		}
		return vals;
	}

	/** Get row dimension.
   @return     m, the number of rows.
	 */

	public int getRowDimension () {
		return m;
	}

	/** Get column dimension.
   @return     n, the number of columns.
	 */

	public int getColumnDimension () {
		return n;
	}

	/** Get a single element.
   @param i    Row index.
   @param j    Column index.
   @return     A(i,j)
   @exception  ArrayIndexOutOfBoundsException
	 */

	public double get (int i, int j) {
		return A[i][j];
	}

	/** Get a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param j0   Initial column index
   @param j1   Final column index
   @return     A(i0:i1,j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */

	public Matrix getMatrix (int i0, int i1, int j0, int j1) {
		Matrix X = new Matrix(i1-i0+1,j1-j0+1);
		double[][] B = X.getArray();
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = j0; j <= j1; j++) {
					B[i-i0][j-j0] = A[i][j];
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
		return X;
	}

	/** Get a submatrix.
   @param r    Array of row indices.
   @param c    Array of column indices.
   @return     A(r(:),c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */

	public Matrix getMatrix (int[] r, int[] c) {
		Matrix X = new Matrix(r.length,c.length);
		double[][] B = X.getArray();
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = 0; j < c.length; j++) {
					B[i][j] = A[r[i]][c[j]];
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
		return X;
	}

	/** Get a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param c    Array of column indices.
   @return     A(i0:i1,c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */

	public Matrix getMatrix (int i0, int i1, int[] c) {
		Matrix X = new Matrix(i1-i0+1,c.length);
		double[][] B = X.getArray();
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = 0; j < c.length; j++) {
					B[i-i0][j] = A[i][c[j]];
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
		return X;
	}

	/** Get a submatrix.
   @param r    Array of row indices.
   @param i0   Initial column index
   @param i1   Final column index
   @return     A(r(:),j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */

	public Matrix getMatrix (int[] r, int j0, int j1) {
		Matrix X = new Matrix(r.length,j1-j0+1);
		double[][] B = X.getArray();
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = j0; j <= j1; j++) {
					B[i][j-j0] = A[r[i]][j];
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
		return X;
	}

	/** Set a single element.
   @param i    Row index.
   @param j    Column index.
   @param s    A(i,j).
   @exception  ArrayIndexOutOfBoundsException
	 */

	public void set (int i, int j, double s) {
		A[i][j] = s;
	}

	/** Set a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param j0   Initial column index
   @param j1   Final column index
   @param X    A(i0:i1,j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */

	public void setMatrix (int i0, int i1, int j0, int j1, Matrix X) {
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = j0; j <= j1; j++) {
					A[i][j] = X.get(i-i0,j-j0);
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/** Set a submatrix.
   @param r    Array of row indices.
   @param c    Array of column indices.
   @param X    A(r(:),c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */

	public void setMatrix (int[] r, int[] c, Matrix X) {
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = 0; j < c.length; j++) {
					A[r[i]][c[j]] = X.get(i,j);
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/** Set a submatrix.
   @param r    Array of row indices.
   @param j0   Initial column index
   @param j1   Final column index
   @param X    A(r(:),j0:j1)
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */

	public void setMatrix (int[] r, int j0, int j1, Matrix X) {
		try {
			for (int i = 0; i < r.length; i++) {
				for (int j = j0; j <= j1; j++) {
					A[r[i]][j] = X.get(i,j-j0);
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/** Set a submatrix.
   @param i0   Initial row index
   @param i1   Final row index
   @param c    Array of column indices.
   @param X    A(i0:i1,c(:))
   @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */

	public void setMatrix (int i0, int i1, int[] c, Matrix X) {
		try {
			for (int i = i0; i <= i1; i++) {
				for (int j = 0; j < c.length; j++) {
					A[i][c[j]] = X.get(i-i0,j);
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException("Submatrix indices");
		}
	}

	/** Matrix transpose.
   @return    A'
	 */

	public Matrix transpose () {
		Matrix X = new Matrix(n,m);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[j][i] = A[i][j];
			}
		}
		return X;
	}

	/** One norm
   @return    maximum column sum.
	 */

	public double norm1 () {
		double f = 0;
		for (int j = 0; j < n; j++) {
			double s = 0;
			for (int i = 0; i < m; i++) {
				s += Math.abs(A[i][j]);
			}
			f = Math.max(f,s);
		}
		return f;
	}



	/** Infinity norm
   @return    maximum row sum.
	 */

	public double normInf () {
		double f = 0;
		for (int i = 0; i < m; i++) {
			double s = 0;
			for (int j = 0; j < n; j++) {
				s += Math.abs(A[i][j]);
			}
			f = Math.max(f,s);
		}
		return f;
	}

	/** Frobenius norm
   @return    sqrt of sum of squares of all elements.
	 */

	public double normF () {
		double f = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				f = Utility.hypot(f,A[i][j]);
			}
		}
		return f;
	}

	/**  Unary minus
   @return    -A
	 */

	public Matrix uminus () {
		Matrix X = new Matrix(m,n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = -A[i][j];
			}
		}
		return X;
	}

	/** C = A + B
   @param B    another matrix
   @return     A + B
	 */

	public Matrix plus (Matrix B) {
		checkMatrixDimensions(B);
		Matrix X = new Matrix(m,n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] + B.A[i][j];
			}
		}
		return X;
	}

	/** A = A + B
   @param B    another matrix
   @return     A + B
	 */

	public Matrix plusEquals (Matrix B) {
		checkMatrixDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = A[i][j] + B.A[i][j];
			}
		}
		return this;
	}

	/** C = A - B
   @param B    another matrix
   @return     A - B
	 */

	public Matrix minus (Matrix B) {
		checkMatrixDimensions(B);
		Matrix X = new Matrix(m,n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] - B.A[i][j];
			}
		}
		return X;
	}

	/** A = A - B
   @param B    another matrix
   @return     A - B
	 */

	public Matrix minusEquals (Matrix B) {
		checkMatrixDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = A[i][j] - B.A[i][j];
			}
		}
		return this;
	}

	/** Element-by-element multiplication, C = A.*B
   @param B    another matrix
   @return     A.*B
	 */

	public Matrix arrayTimes (Matrix B) {
		checkMatrixDimensions(B);
		Matrix X = new Matrix(m,n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] * B.A[i][j];
			}
		}
		return X;
	}

	/** Element-by-element multiplication in place, A = A.*B
   @param B    another matrix
   @return     A.*B
	 */

	public Matrix arrayTimesEquals (Matrix B) {
		checkMatrixDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = A[i][j] * B.A[i][j];
			}
		}
		return this;
	}

	/** Element-by-element right division, C = A./B
   @param B    another matrix
   @return     A./B
	 */

	public Matrix arrayRightDivide (Matrix B) {
		checkMatrixDimensions(B);
		Matrix X = new Matrix(m,n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] / B.A[i][j];
			}
		}
		return X;
	}

	/** Element-by-element right division in place, A = A./B
   @param B    another matrix
   @return     A./B
	 */

	public Matrix arrayRightDivideEquals (Matrix B) {
		checkMatrixDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = A[i][j] / B.A[i][j];
			}
		}
		return this;
	}

	/** Element-by-element left division, C = A.\B
   @param B    another matrix
   @return     A.\B
	 */

	public Matrix arrayLeftDivide (Matrix B) {
		checkMatrixDimensions(B);
		Matrix X = new Matrix(m,n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = B.A[i][j] / A[i][j];
			}
		}
		return X;
	}

	/** Element-by-element left division in place, A = A.\B
   @param B    another matrix
   @return     A.\B
	 */

	public Matrix arrayLeftDivideEquals (Matrix B) {
		checkMatrixDimensions(B);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = B.A[i][j] / A[i][j];
			}
		}
		return this;
	}

	/** Multiply a matrix by a scalar, C = s*A
   @param s    scalar
   @return     s*A
	 */

	public Matrix times (double s) {
		Matrix X = new Matrix(m,n);
		double[][] C = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = s*A[i][j];
			}
		}
		return X;
	}

	/** Multiply a matrix by a scalar in place, A = s*A
   @param s    scalar
   @return     replace A by s*A
	 */

	public Matrix timesEquals (double s) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = s*A[i][j];
			}
		}
		return this;
	}

	/** Linear algebraic matrix multiplication, A * B
   @param B    another matrix
   @return     Matrix product, A * B
   @exception  IllegalArgumentException Matrix inner dimensions must agree.
	 */

	public Matrix times (Matrix B) {
		if (B.m != n) {
			throw new IllegalArgumentException("Matrix inner dimensions must agree.");
		}
		Matrix X = new Matrix(m,B.n);
		double[][] C = X.getArray();
		double[] Bcolj = new double[n];
		for (int j = 0; j < B.n; j++) {
			for (int k = 0; k < n; k++) {
				Bcolj[k] = B.A[k][j];
			}
			for (int i = 0; i < m; i++) {
				double[] Arowi = A[i];
				double s = 0;
				for (int k = 0; k < n; k++) {
					s += Arowi[k]*Bcolj[k];
				}
				C[i][j] = s;
			}
		}
		return X;
	}





	/** Matrix trace.
   @return     sum of the diagonal elements.
	 */

	public double trace () {
		double t = 0;
		for (int i = 0; i < Math.min(m,n); i++) {
			t += A[i][i];
		}
		return t;
	}

	/** Generate matrix with random elements
   @param m    Number of rows.
   @param n    Number of colums.
   @return     An m-by-n matrix with uniformly distributed random elements.
	 */

	public static Matrix random (int m, int n) {
		Matrix A = new Matrix(m,n);
		double[][] X = A.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				X[i][j] = Math.random();
			}
		}
		return A;
	}

	/** Generate identity matrix
   @param m    Number of rows.
   @param n    Number of colums.
   @return     An m-by-n matrix with ones on the diagonal and zeros elsewhere.
	 */

	public static Matrix identity (int m, int n) {
		Matrix A = new Matrix(m,n);
		double[][] X = A.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				X[i][j] = (i == j ? 1.0 : 0.0);
			}
		}
		return A;
	}


	/** Print the matrix to stdout.   Line the elements up in columns
	 * with a Fortran-like 'Fw.d' style format.
   @param w    Column width.
   @param d    Number of digits after the decimal.
	 */

	public void print (int w, int d) {
		print(new PrintWriter(System.out,true),w,d); }

	/** Print the matrix to the output stream.   Line the elements up in
	 * columns with a Fortran-like 'Fw.d' style format.
   @param output Output stream.
   @param w      Column width.
   @param d      Number of digits after the decimal.
	 */

	public void print (PrintWriter output, int w, int d) {
		DecimalFormat format = new DecimalFormat();
		format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
		format.setMinimumIntegerDigits(1);
		format.setMaximumFractionDigits(d);
		format.setMinimumFractionDigits(d);
		format.setGroupingUsed(false);
		print(output,format,w+2);
	}

	/** Print the matrix to stdout.  Line the elements up in columns.
	 * Use the format object, and right justify within columns of width
	 * characters.
	 * Note that is the matrix is to be read back in, you probably will want
	 * to use a NumberFormat that is set to US Locale.
   @param format A  Formatting object for individual elements.
   @param width     Field width for each column.
   @see java.text.DecimalFormat#setDecimalFormatSymbols
	 */

	public void print (NumberFormat format, int width) {
		print(new PrintWriter(System.out,true),format,width); }

	// DecimalFormat is a little disappointing coming from Fortran or C's printf.
	// Since it doesn't pad on the left, the elements will come out different
	// widths.  Consequently, we'll pass the desired column width in as an
	// argument and do the extra padding ourselves.

	/** Print the matrix to the output stream.  Line the elements up in columns.
	 * Use the format object, and right justify within columns of width
	 * characters.
	 * Note that is the matrix is to be read back in, you probably will want
	 * to use a NumberFormat that is set to US Locale.
   @param output the output stream.
   @param format A formatting object to format the matrix elements 
   @param width  Column width.
   @see java.text.DecimalFormat#setDecimalFormatSymbols
	 */

	public void print (PrintWriter output, NumberFormat format, int width) {
		output.println();  // start on new line.
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				String s = format.format(A[i][j]); // format the number
				int padding = Math.max(1,width-s.length()); // At _least_ 1 space
				for (int k = 0; k < padding; k++)
					output.print(' ');
				output.print(s);
			}
			output.println();
		}
		output.println();   // end with blank line.
	}

	/** Read a matrix from a stream.  The format is the same the print method,
	 * so printed matrices can be read back in (provided they were printed using
	 * US Locale).  Elements are separated by
	 * whitespace, all the elements for each row appear on a single line,
	 * the last row is followed by a blank line.
   @param input the input stream.
	 */

	public static Matrix read (BufferedReader input, int header) throws java.io.IOException {
		StringTokenizer ST;
		double[][] A ;
		String str;
		int n=0,m=0;
		if ((str=input.readLine())!=null){
			ST = new StringTokenizer(str,",");
			n=(int)Double.parseDouble(ST.nextToken());
			m=(int)Double.parseDouble(ST.nextToken());
		}
		A = new double [n][m];
		for(int i=0; i<n;++i){
			int j=0;
			str=input.readLine();
			ST = new StringTokenizer(str,",");
			while(ST.hasMoreTokens()){
				A[i][j] = Double.parseDouble(ST.nextToken());
				j++;
			}
		}

	
		return new Matrix(A);
	}
	
	public static Matrix readWSpace (BufferedReader input, int header) throws java.io.IOException {
		StringTokenizer ST;
		double[][] A ;
		String str;
		int n=0,m=0;
		if ((str=input.readLine())!=null){
			ST = new StringTokenizer(str," ");
			n=(int)Double.parseDouble(ST.nextToken());
			m=(int)Double.parseDouble(ST.nextToken());
		}
		A = new double [n][m];
		for(int i=0; i<n;++i){
			int j=0;
			str=input.readLine();
			ST = new StringTokenizer(str," ");
			while(ST.hasMoreTokens()){
				A[i][j] = Double.parseDouble(ST.nextToken());
				j++;
			}
		}

	
		return new Matrix(A);
	}
	public void writeMatrixIntoFile(String fileName,int header,int integ){
		PrintWriter out=null;
		DecimalFormat df = new DecimalFormat("#.####");
		try {
			out = new PrintWriter(new FileWriter(fileName));
			if (header==1){
				out.println(this.getRowDimension()+","+this.getColumnDimension());
			}
			for(int i=0; i<this.getRowDimension();++i){
				for(int j=0; j<this.getColumnDimension();++j){
					if(integ==1)
						out.print((int)A[i][j]+",");
					else if(integ==2)
						out.print((byte)A[i][j]+",");
					else
						out.print(df.format(A[i][j])+",");
				}
				out.println();	
			}

			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public void writeMatrixIntoByteFile(String fileName,int header){
		OutputStream out=null;

		try {
			out = new FileOutputStream(fileName);
			if (header==1){
				byte[] a ={(byte)this.getRowDimension(), (byte)this.getColumnDimension()};

				out.write(a);
			}
			for(int i=0; i<this.getRowDimension();++i){
				for(int j=0; j<this.getColumnDimension();++j){
					byte[] a ={(byte)A[i][j]};

					out.write(a);


				}
			}

			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}


	public int maxMatrixEntry(){
		int maxEntry;
		int i, j;

		maxEntry = (int)A[0][0];
		for (i = 0; i < m; i++)
			for (j = 0; j < n; j++)
				if (maxEntry < (int)A[i][j])
					maxEntry = (int)A[i][j];
		return maxEntry;
	}

	public void relabelComponents(){
		int maxLabel = this.maxMatrixEntry();
		int currLabel;
		int relabels[] = new int [maxLabel+1]; 
		int i, j;

		for (i = 0; i <= maxLabel; i++)
			relabels[i] = -1;
		currLabel = 1;
		for (i = 0 ; i < this.getRowDimension(); i++)
			for (j = 0; j < this.getColumnDimension(); j++)
				if (this.get(i,j) > 0 ) 
					if (relabels[(int)this.get(i, j)] == -1)
						relabels[ (int)this.get(i, j)] = currLabel++;
		for (i = 0 ; i <  this.getRowDimension(); i++)
			for (j = 0; j < this.getColumnDimension(); j++)
				if ( this.get(i,j) > 0)
					this.set(i, j, relabels[ (int) this.get(i, j) ]);
	}


	public static Matrix generateCircularWindow(int winSize){
		Matrix M = new Matrix(winSize,winSize,1);
		int i, j, space;

		space = (winSize / 4);
		for (i = 0; i < space; i++)
			for (j = 0; j < space - i; j++){
				M.set(i, j, 0);
				M.set(M.getRowDimension() - 1 - i, j,0);
				M.set(i, M.getColumnDimension()-1-j,0);
				M.set(M.getRowDimension() - 1 - i, M.getColumnDimension()-1-j,0);
			}
		return M;
	}
	public Matrix markRegionOfInterest(int regionOfInterest) {
		// TODO Auto-generated method stub
		Matrix result = new Matrix(this.m,this.n,Utility.BACKGROUND);
		int i, j;

		for (i = 0; i < this.m; i++)
			for (j = 0; j < this.n; j++)
				if (this.A[i][j] == regionOfInterest)
					result.A[i][j] = 1;
		return result;
	}

	public Matrix fourConnectivity() {
		// TODO Auto-generated method stub
		Matrix visited = new Matrix(this.m,this.n,Utility.UNVISITED);
		Matrix regions = new Matrix(this.m,this.n,Utility.BACKGROUND);
		int queueX[] = new int[this.m*this.n]; 
		int queueY[] = new int[this.m*this.n]; 
		int i, j, x, y, queueStart, queueEnd;
		int label = 0;


		for (i = 0; i < this.m; i++)
			for (j = 0; j < this.n; j++){
				if (this.A[i][j] == Utility.BACKGROUND)
					continue;
				if (visited.A[i][j] == Utility.VISITED)
					continue;

				label++;
				queueX[0] = i;		queueY[0] = j;
				queueStart = 0;		queueEnd = 1;
				visited.A[i][j] = Utility.VISITED;

				while (queueEnd > queueStart){
					x = queueX[queueStart];		y = queueY[queueStart];
					regions.A[x][y] = label;
					if (x-1 >= 0 && visited.A[x-1][y] != Utility.VISITED && this.A[x-1][y] != Utility.BACKGROUND){
						visited.A[x-1][y] = Utility.VISITED;
						queueX[queueEnd] = x-1;		queueY[queueEnd] = y;		queueEnd++;
					}
					if (x+1 < this.m && visited.A[x+1][y] != Utility.VISITED && this.A[x+1][y] != Utility.BACKGROUND){
						visited.A[x+1][y] = Utility.VISITED;
						queueX[queueEnd] = x+1;		queueY[queueEnd] = y;		queueEnd++;
					}
					if (y-1 >= 0 && visited.A[x][y-1] != Utility.VISITED && this.A[x][y-1] != Utility.BACKGROUND){
						visited.A[x][y-1] = Utility.VISITED;
						queueX[queueEnd] = x;		queueY[queueEnd] = y-1;		queueEnd++;
					}
					if (y+1 < this.n && visited.A[x][y+1] != Utility.VISITED && this.A[x][y+1] != Utility.BACKGROUND){
						visited.A[x][y+1] = Utility.VISITED;
						queueX[queueEnd] = x;		queueY[queueEnd] = y+1;		queueEnd++;
					}
					queueStart++;
				}
			}
		this.setUnsegmentedRegionNo(label);
		regions.setUnsegmentedRegionNo(label);
		return regions;
	}

	public Matrix eightConnectivity() {
		// TODO Auto-generated method stub
		Matrix visited = new Matrix(this.m,this.n,Utility.UNVISITED);
		Matrix regions = new Matrix(this.m,this.n,Utility.BACKGROUND);
		int queueX[] = new int[this.m*this.n]; 
		int queueY[] = new int[this.m*this.n]; 
		int i, j, x, y, queueStart, queueEnd;
		int label = 0;


		for (i = 0; i < this.m; i++)
			for (j = 0; j < this.n; j++){
				if (this.A[i][j] == Utility.BACKGROUND)
					continue;
				if (visited.A[i][j] == Utility.VISITED)
					continue;

				label++;
				queueX[0] = i;		queueY[0] = j;
				queueStart = 0;		queueEnd = 1;
				visited.A[i][j] = Utility.VISITED;

				while (queueEnd > queueStart){
					x = queueX[queueStart];		y = queueY[queueStart];
					regions.A[x][y] = label;

					if (x-1 >= 0 && visited.A[x-1][y] != Utility.VISITED && this.A[x-1][y] != Utility.BACKGROUND){
						visited.A[x-1][y] = Utility.VISITED;
						queueX[queueEnd] = x-1;		queueY[queueEnd] = y;		queueEnd++;
					}
					if (x+1 < this.m && visited.A[x+1][y] != Utility.VISITED && this.A[x+1][y] != Utility.BACKGROUND){
						visited.A[x+1][y] = Utility.VISITED;
						queueX[queueEnd] = x+1;		queueY[queueEnd] = y;		queueEnd++;
					}
					if (y-1 >= 0 && visited.A[x][y-1] != Utility.VISITED && this.A[x][y-1] != Utility.BACKGROUND){
						visited.A[x][y-1] = Utility.VISITED;
						queueX[queueEnd] = x;		queueY[queueEnd] = y-1;		queueEnd++;
					}
					if (y+1 < this.n && visited.A[x][y+1] != Utility.VISITED && this.A[x][y+1] != Utility.BACKGROUND){
						visited.A[x][y+1] = Utility.VISITED;
						queueX[queueEnd] = x;		queueY[queueEnd] = y+1;		queueEnd++;
					}
					if (x-1 >= 0 && y-1 >=0 && visited.A[x-1][y-1] != Utility.VISITED && this.A[x-1][y-1] != Utility.BACKGROUND){
						visited.A[x-1][y-1] = Utility.VISITED;
						queueX[queueEnd] = x-1;		queueY[queueEnd] = y-1;		queueEnd++;
					}
					if (x-1 >= 0 && y+1 < this.n && visited.A[x-1][y+1] != Utility.VISITED && this.A[x-1][y+1] != Utility.BACKGROUND){
						visited.A[x-1][y+1] = Utility.VISITED;
						queueX[queueEnd] = x-1;		queueY[queueEnd] = y+1;		queueEnd++;
					}
					if (x+1 < this.m && y-1 >=0 && visited.A[x+1][y-1] != Utility.VISITED && this.A[x+1][y-1] != Utility.BACKGROUND){
						visited.A[x+1][y-1] = Utility.VISITED;
						queueX[queueEnd] = x+1;		queueY[queueEnd] = y-1;		queueEnd++;
					}
					if (x+1 < this.m && y+1 < this.n && visited.A[x+1][y+1] != Utility.VISITED && this.A[x+1][y+1] != Utility.BACKGROUND){
						visited.A[x+1][y+1] = Utility.VISITED;
						queueX[queueEnd] = x+1;		queueY[queueEnd] = y+1;		queueEnd++;
					}
					queueStart++;
				}
			}
		this.setUnsegmentedRegionNo(label);
		regions.setUnsegmentedRegionNo(label);
		return regions;
	}

	public int countMatrixOccurrences( int value){
		int occ = 0, i, j;

		for (i = 0; i < this.m; i++)
			for (j = 0; j < this.n; j++)
				if (this.A[i][j] == value)
					occ++;
		return occ;
	}

	public void findAdjacents(Matrix prev, Matrix curr, int connectivity) {
		// TODO Auto-generated method stub		
		int i, j;

		for (i = 0; i < curr.getRowDimension(); i++)
			for (j = 0; j < curr.getColumnDimension(); j++)
				if (curr.get(i, j) > 0){	// if this is a segmented pixel
					if ((i - 1 >= 0) && (prev.get(i-1, j) > 0))
						this.A[ (int)curr.get(i, j) ][ (int)prev.get(i-1, j) ] = 1;
					if ((i + 1 < curr.getRowDimension()) && (prev.get(i+1, j) > 0))
						this.A[ (int)curr.get(i, j) ][ (int)prev.get(i+1, j) ] = 1;
					if ((j - 1 >= 0) && ( prev.get(i, j-1) > 0))
						this.A[ (int)curr.get(i, j) ][ (int)prev.get(i, j-1) ] = 1;
					if ((j + 1 < curr.getColumnDimension()) && ((int)prev.get(i, j+1) > 0))
						this.A[ (int)curr.get(i, j) ][ (int)prev.get(i, j+1) ] = 1;
					if (connectivity == Utility.EIGHT){
						if ((i - 1 >= 0) && (j - 1 >= 0) && ((int)prev.get(i-1, j-1) > 0))
							this.A[ (int)curr.get(i, j) ][ (int)prev.get(i-1, j-1) ] = 1;
						if ((i - 1 >= 0) && (j + 1 < curr.getColumnDimension()) && ((int)prev.get(i-1, j+1) > 0))
							this.A[ (int)curr.get(i, j) ][ (int)prev.get(i-1, j+1) ] = 1;
						if ((i + 1 < curr.getRowDimension()) && (j - 1 >= 0) && ((int)prev.get(i+1, j-1) > 0))
							this.A[ (int)curr.get(i, j) ][ (int)prev.get(i+1, j-1) ] = 1;
						if ((i + 1 < curr.getRowDimension()) && (j + 1 < curr.getColumnDimension()) && ((int)prev.get(i+1, j+1) > 0))
							this.A[ (int)curr.get(i, j) ][ (int)prev.get(i+1, j+1) ] = 1;
					}
				}
	}




	public void removeSmallNumberOfObjects(int areaThr, Matrix objects) {
		int rno = this.maxMatrixEntry();
		int i, j;
		int noOfObjects[] = new int[rno+1]; 
		int visited[] = new int[objects.maxMatrixEntry()+1];
		for (i = 0; i < this.getRowDimension(); i++)
			for (j = 0; j < this.getColumnDimension(); j++){
				int ono = (int)objects.get(i, j);
				if((int)this.get(i, j)>0 && visited[ono]!=1){
					noOfObjects[ (int)this.get(i, j) ]++;
					visited[ono]=1;
				}
			}
		for (i = 0; i < this.getRowDimension(); i++)
			for (j = 0; j < this.getColumnDimension(); j++)
				if (noOfObjects[ (int)this.get(i, j) ] < areaThr)
					A[i][j] = Utility.UNSEGMENTED;	


	}


	public void removeHoles(int connectivity) {
		Matrix temp = this.markRegionOfInterest(Utility.UNSEGMENTED);
		Matrix unsegmentedRegions=null;
		int neigh[], i, j;

		if (connectivity == Utility.FOUR)
			unsegmentedRegions = temp.fourConnectivity();
		else if (connectivity == Utility.EIGHT)
			unsegmentedRegions = temp.eightConnectivity();
		else{
			System.out.println("\nError: Connectivity could be either four or eight\n\n");
			System.exit(1);
		}
		neigh = new int [unsegmentedRegions.getUnsegmentedRegionNo()+1]; // initialized with zero

		for (i = 0; i < unsegmentedRegions.getRowDimension(); i++)
			for (j = 0; j < unsegmentedRegions.getColumnDimension(); j++){
				if (neigh [ (int)unsegmentedRegions.get(i,j) ] == -1)
					continue;
				if ((i <= 0) || (j <= 0) || 
						(i >= unsegmentedRegions.getRowDimension() - 1) || (j >= unsegmentedRegions.getColumnDimension() - 1))
					continue;
				if (this.get(i-1, j) != Utility.UNSEGMENTED){
					if (neigh[ (int)unsegmentedRegions.get(i,j) ] == 0)
						neigh[ (int)unsegmentedRegions.get(i,j) ] = (int) A[i-1][j];
					else if (neigh[ (int)unsegmentedRegions.get(i,j)] != A[i-1][j])
						neigh[ (int)unsegmentedRegions.get(i,j)] = -1;
				}
				if (this.get(i+1, j) != Utility.UNSEGMENTED){
					if (neigh[ (int)unsegmentedRegions.get(i,j) ] == 0)
						neigh[ (int)unsegmentedRegions.get(i,j) ] = (int) A[i+1][j];
					else if (neigh[ (int)unsegmentedRegions.get(i,j)] != A[i+1][j])
						neigh[ (int)unsegmentedRegions.get(i,j)] = -1;
				}
				if (this.get(i, j-1) != Utility.UNSEGMENTED){
					if (neigh[ (int)unsegmentedRegions.get(i,j) ] == 0)
						neigh[ (int)unsegmentedRegions.get(i,j) ] = (int)A[i][j-1];
					else if (neigh[ (int)unsegmentedRegions.get(i,j) ] != A[i][j-1])
						neigh[ (int)unsegmentedRegions.get(i,j) ] = -1;
				}
				if (this.get(i, j+1) != Utility.UNSEGMENTED){
					if (neigh[ (int)unsegmentedRegions.get(i,j) ] == 0)
						neigh[ (int)unsegmentedRegions.get(i,j) ] = (int)A[i][j+1];
					else if (neigh[ (int)unsegmentedRegions.get(i,j) ] != A[i][j+1])
						neigh[ (int)unsegmentedRegions.get(i,j) ] = -1;
				}
			}
		for (i = 0; i < unsegmentedRegions.getRowDimension(); i++)
			for (j = 0; j < unsegmentedRegions.getColumnDimension(); j++)
				if (neigh[ (int)unsegmentedRegions.get(i,j) ] > 0)
					A[i][j] = neigh[ (int)unsegmentedRegions.get(i,j) ];

	}



	Matrix findAllConnectedComponents(int K, int connectivity){
		Matrix finalRegions, currentMap;
		Matrix regions[] = new Matrix[K]; //(MATRIX *) malloc(K * sizeof(MATRIX));
		int labels[] = new int[K]; //(int *) malloc(K * sizeof(int));
		int currLabel, i, j, k;

		for (k = 0; k < K; k++){
			currentMap = this.markRegionOfInterest(k);
			if (connectivity == Utility.FOUR){
				regions[k] = currentMap.fourConnectivity();
				labels[k]=regions[k].getUnsegmentedRegionNo();
			}
			else if (connectivity == Utility.EIGHT){
				regions[k] = currentMap.eightConnectivity();
				labels[k]=regions[k].getUnsegmentedRegionNo();
			}
			else{

			}
		}


		finalRegions = new Matrix(this.getRowDimension(),this.getColumnDimension(),Utility.BACKGROUND);
		currLabel = 0;
		for (k = 0; k < K; k++){
			for (i = 0; i < finalRegions.getRowDimension(); i++)
				for (j = 0; j < finalRegions.getColumnDimension(); j++){
					if (regions[k].get(i, j) == Utility.BACKGROUND)
						continue;
					finalRegions.set(i, j, currLabel + regions[k].get(i, j));
				}
			currLabel += labels[k];
		}

		return finalRegions;
	}






	public Matrix normalize() {
		// TODO Auto-generated method stub


		int p,q,i,j;
		double length;
		Matrix result = this.copy();
		for(p=0; p<result.getColumnDimension(); p++) // Loops over all vectors
		{
			// Calculates the normalization factor
			length = 0;
			for(i=0; i<result.getRowDimension(); i++)
				for(j=0; j<result.getRowDimension(); j++)
					length += result.get(i, p)*result.get(j, p)*this.get(i, j);

			length = Math.sqrt(length);

			// Normalizes the vector
			if (length!=0d)
				for(q=0; q<result.getRowDimension(); q++)
					result.set(q, p, result.get(q, p)/length);
			else
				System.out.println("Warning(orthonormalize):"+(p+1)+". Vector has length null");
		}
		return result;


	}





	public void enlargeByColumn(int numOfColumns){
		double[][] newData = new double[this.getRowDimension()][this.getColumnDimension()+numOfColumns];

		for(int i=0; i<this.getRowDimension(); i++){
			for(int j=0; j<this.getColumnDimension(); j++)
				newData[i][j] = this.A[i][j];
			for(int j=this.getColumnDimension(); j<this.getColumnDimension()+numOfColumns; j++)
				newData[i][j] = 0;
		}

		this.n = this.getColumnDimension() + numOfColumns;
		this.A = newData;
	}



	public static Matrix deleteTheColumnsWithAllZeros(Matrix M){
		int beginCol=1;


		int endCol=M.getColumnDimension()-1;
		boolean advancePrevCol = true;
		while(advancePrevCol){
			for(int i=0; i<M.getRowDimension(); i++)
				if((int)M.get(i, endCol) != 0){
					advancePrevCol = false;
					break;
				}

			if(advancePrevCol)
				endCol--;
		}

		Matrix ret = new Matrix(M.getRowDimension(), endCol-beginCol+1);
		for(int row=0; row<M.getRowDimension(); row++){
			for(int i=0, j=beginCol; j<=endCol; i++,j++){
				ret.set(row, i, M.get(row, j));
			}
		}

		return ret;
	}


	private void checkMatrixDimensions (Matrix B) {
		if (B.m != m || B.n != n) {
			throw new IllegalArgumentException("Matrix dimensions must agree.");
		}
	}

	public int getUnsegmentedRegionNo() {
		return unsegmentedRegionNo;
	}

	public void setUnsegmentedRegionNo(int unsegmentedRegionNo) {
		this.unsegmentedRegionNo = unsegmentedRegionNo;
	}



	public Matrix form2DMatrix(int row, int column){
		Matrix newM;
		int i, j, c;

		if (row * column != this.getRowDimension()){
			System.out.println("\nError: dimensions mismatch\n\n");
			System.exit(1);
		}
		newM = new Matrix(row,column);
		c = 0;
		for (i = 0; i < row; i++)
			for (j = 0; j < column; j++)
				newM.A[i][j] = this.A[c++][0];
		return newM;
	}

	public double maxAbsoluteDifferenceBetweenMatrices(Matrix B){
		int i, j;
		double maxAbsEntry = Math.abs(this.A[0][0] - B.A[0][0]);

		for (i = 0; i < this.getRowDimension(); i++)
			for (j = 0; j < this.getColumnDimension(); j++)
				if (maxAbsEntry < Math.abs(this.A[i][j] - B.A[i][j]))
					maxAbsEntry = Math.abs(this.A[i][j] - B.A[i][j]);
		return maxAbsEntry;
	}

	public void initializeMatrix (double x) {

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = x;
			}
		}

	}

	public Matrix computeCorrelationMatrix(){
		int i, j;
		Matrix S = this.computeCovarianceMatrix();
		Matrix R = new Matrix (S.getRowDimension(),S.getColumnDimension());
		for (i = 0; i < S.getRowDimension(); i++)
			for (j = 0; j < S.getColumnDimension(); j++)
				if (Math.abs(Math.sqrt(S.get(i, i)) * Math.sqrt(S.get(j, j))) > Utility.ZERO)
					R.set(i, j, S.get(i, j)/ (Math.sqrt(S.get(i, i)) * Math.sqrt(S.get(j, j)))) ;
				else
					R.set(i, j,0);

		return R;
	}

	public Matrix computeCovarianceMatrix(){
		int i, j, t;
		Matrix M = this.computeMeanMatrix();
		Matrix S = new Matrix(this.getColumnDimension(),this.getColumnDimension(),0.0);


		for (i = 0; i < S.getRowDimension(); i++)
			for (j = 0; j < S.getColumnDimension(); j++)
				for (t = 0; t < this.getRowDimension(); t++)
					S.set(i, j, S.get(i, j)+ (this.get(t, i) - M.get(0, i)) * (this.get(t, j) - M.get(0, j)) );
		for (i = 0; i < S.getRowDimension(); i++)
			for (j = 0; j < S.getColumnDimension(); j++)
				S.set(i, j, S.get(i, j)/ (this.getRowDimension() - 1));

		return S;
	}

	public Matrix computeMeanMatrix(){
		int i, j;
		Matrix M = new Matrix(1,this.getColumnDimension(),0.0);


		for (i = 0; i < this.getRowDimension(); i++)
			for (j = 0; j < this.getColumnDimension(); j++)
				M.set(0, j, M.get(0, j)+this.get(i,j));

		for (j = 0; j < M.getColumnDimension(); j++)
			M.set(0, j, M.get(0, j)/this.getRowDimension());
		return M;
	}
	// for symmetric matrices
	public Matrix computeEigenValues(Matrix V){
		int i, j;

		Matrix D;

		if (this.getRowDimension() != this.getColumnDimension()){
			System.out.println("\nError: Matrix should be square for eigenvalue computation");
			System.exit(1);
		}
		for (i = 0; i < this.getRowDimension(); i++)
			for (j = i+1; j < this.getColumnDimension(); j++)
				if (this.get(i, j) != this.get(j, i)){
					System.out.println("\nError: Matrix should be symmetric for eigenvalue computation");
					System.exit(1);
				}
		D = new Matrix(1,this.getRowDimension());
		int [] nrot = new int[1];
		Utility.jacobi(this.A,this.getRowDimension(),D.A[0],V.A,nrot);
		Utility.eigsrt(D.A[0],V.A,D.getColumnDimension());
		return D;
	}

	public Matrix multiplyMatrix( Matrix B){
		Matrix result;
		int i, j, k;

		if (this.getColumnDimension() != B.getRowDimension()){
			System.out.println("\nError: Matrix dimensions do not match in matrix multiplication");
			System.exit(1);
		}
		result = new Matrix(this.getRowDimension(),B.getColumnDimension());
		for (i = 0; i < this.getRowDimension(); i++)
			for (j = 0; j < B.getColumnDimension(); j++){
				result.A[i][j] = 0.0;
				for (k = 0; k < this.getColumnDimension(); k++)
					result.A[i][j] += this.A[i][k] * B.A[k][j];
			}
		return result;
	}

	public double maxMatrixEntry(int whichColumn){
		double maxEntry;
		int i, j;

		if (whichColumn == -1){
			maxEntry = this.A[0][0];
			for (i = 0; i < this.getRowDimension(); i++)
				for (j = 0; j < this.getColumnDimension(); j++)
					if (maxEntry < this.A[i][j])
						maxEntry = this.A[i][j];
			return maxEntry;
		}
		maxEntry = this.A[0][whichColumn];
		for (i = 0; i < this.getRowDimension(); i++)
			if (maxEntry < this.A[i][whichColumn])
				maxEntry = this.A[i][whichColumn];
		return maxEntry;
	}
	public double minMatrixEntry( int whichColumn){
		double minEntry;
		int i, j;

		if (whichColumn == -1){
			minEntry = this.A[0][0];
			for (i = 0; i < this.getRowDimension(); i++)
				for (j = 0; j < this.getColumnDimension(); j++)
					if (minEntry > this.A[i][j])
						minEntry = this.A[i][j];
			return minEntry;
		}
		minEntry = this.A[0][whichColumn];
		for (i = 0; i < this.getRowDimension(); i++)
			if (minEntry > this.A[i][whichColumn])
				minEntry = this.A[i][whichColumn];
		return minEntry;
	}

	public Matrix inverseMatrix(){
		Matrix inv, temp;
		int i;

		if (this.getRowDimension() != this.getColumnDimension()){
			System.out.println("Error: Matrix should be square for inverse computation");
			System.exit(1);
		}
		temp = new Matrix(this.getRowDimension(), this.getColumnDimension()); 

		temp = this.copy();
		inv = new Matrix(this.getRowDimension(), this.getColumnDimension(),0.0); 

		for (i = 0; i < this.getRowDimension(); i++)
			inv.set(i, i, 1);
		Utility.gaussj(this.A,this.getRowDimension(),inv.A,inv.getColumnDimension());
		this.setMatrix(temp); //copyMatrixD(&M,temp);

		return inv;
	}
	private void setMatrix(Matrix temp) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = temp.get(i, j);
			}
		}

	}

	public double mean() {
		double mean=0;

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				mean = mean + A[i][j];
			}
		}
		return mean/(m*n);
	}

	public double std(double avg) {

		double stdev=0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				stdev += Utility.SQUARE(A[i][j] - avg);
			}
		}
		stdev = Math.sqrt(stdev / ((m*n) - 1));
		return stdev;
	}


	public Matrix findEdges(){
		
		Matrix edge = new Matrix(m,n,0);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if ((i + 1 < m) && ( A[i][j] != A[i+1][j])){
					edge.set(i,j,1);
				}
				if ((j + 1 < n) && (A[i][j] != A[i][j+1])){
					edge.set(i,j,1);
				}
			}
		}
		return edge;
	}


	public void bwdilate2(Matrix se) {
		Matrix R = new Matrix(m,n,0);
		int i, j, k1, k2;
		int msize = se.getRowDimension() / 2;

		if (se.getRowDimension() % 2 == 0){	// if even, exit
			System.out.println("\nError: Size of the structural element should be odd\n\n");
			System.exit(1);
		}

		for (i = 0; i < m; i++)
			for (j = 0; j < n; j++){
				int found = 0;
				if (A[i][j] != 0 && A[i][j] != 1){
					System.out.println("\nError: Only BW images can be dilated\n\n");
					System.exit(1);
				}
				for (k1 = -msize; k1 <= msize; k1++){
					for (k2 = -msize; k2 <= msize; k2++){
						int x = i + k1; 
						int y = j + k2;
						if ((x < 0) || (x >= m) || (y < 0) || (y >= n))
							continue;
						if (se.get(k1 + msize, k2 + msize)>0 && (A[x][y] == 1)){
							R.set(i, j, 1);
							found = 1;
							break;
						}
					}
					if (found==1)
						break;
				}
			}
		for (i = 0; i < m; i++)
			for (j = 0; j < n; j++)
				A[i][j]=R.get(i, j);
	}






}
