 package matrix;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.StringTokenizer;

/**
 * Holds the data for an integer matrix and allows various operations.
 * @author S. Tuncer Erdogan
 */
public class IntMatrix{
	public int[][] data;
	public int row; 
	public int column;
	
	/**
	 * Default constructor, does nothing.
	 */
	public IntMatrix(){
	}
	
	/**
	 * Allocates a data member with given size 
	 * @param row
	 * @param column
	 */
	public IntMatrix(int row, int column){
		init(row,column);
	}
	
	/**
	 * Copies the contents of the matrix to the current object.
	 * If the matrix is not initialized (row=0, column=0), the matrix is initialized
	 * with the same size as the matrix to be coppied.
	 * Otherwise error message is printed and exited.
	 * @param copy IntMatrix to be copied
	 */
	public IntMatrix copyMatrix(IntMatrix copy){
		init(copy.row, copy.column);
		
		for (int i = 0; i < copy.row; i++){
			for (int j = 0; j < copy.column; j++){
				this.data[i][j] = copy.data[i][j];
			}
		}
		
		return this;
	}
	
	public void copyMatrix(int[][] intArray){
		int r = intArray.length;
		int w = intArray[0].length;
		
		init(r, w);
		
		for (int i = 0; i < r; i++){
			for (int j = 0; j < w; j++){
				this.data[i][j] = intArray[i][j];
			}
		}		
	}
	
	public void copyMatrix(IntMatrix copy, int offsetX, int offsetY){
		for(int i=0; i<copy.row; i++){
			for(int j=0; j<copy.column; j++){
				this.data[i+offsetX][j+offsetY] = copy.data[i][j]; 
			}
		}
	}
	
	public static IntMatrix deleteTheColumnsWithAllZeros(IntMatrix A){
		int beginCol=1;
		/*boolean advanceNextCol = true;
		while(advanceNextCol){
			for(int i=0; i<A.row; i++)
				if(A.data[i][beginCol] != 0){
					advanceNextCol = false;
					break;
				}
			
			if(advanceNextCol)
				beginCol++;
		}*/
		
		int endCol=A.column-1;
		boolean advancePrevCol = true;
		while(advancePrevCol){
			for(int i=0; i<A.row; i++)
				if(A.data[i][endCol] != 0){
					advancePrevCol = false;
					break;
				}
			
			if(advancePrevCol)
				endCol--;
		}
		
		IntMatrix ret = new IntMatrix(A.row, endCol-beginCol+1);
		for(int row=0; row<A.row; row++){
			for(int i=0, j=beginCol; j<=endCol; i++,j++){
				ret.data[row][i] = A.data[row][j];
			}
		}
		
		return ret;
	}
	
	/**
	 * Initializes the matrix with the given value, setting all elements of the
	 * matrix to the value.
	 * @param value
	 */
	public void initializeMatrix(int value){
		for(int i=0; i<this.row; i++){
			for(int j=0;j<this.column; j++){
				data[i][j] = value;
			}
		}
	}
	
	public void enlargeByColumn(int numOfColumns){
		int[][] newData = new int[row][column+numOfColumns];
		
		for(int i=0; i<row; i++){
			for(int j=0; j<column; j++)
				newData[i][j] = data[i][j];
			for(int j=column; j<column+numOfColumns; j++)
				newData[i][j] = 0;
		}
		
		this.column = this.column + numOfColumns;
		data = newData;
	}
	
	/**
	 * Returns the specified sub-matrix of the current matrix.
	 * @param offsetRow row number to start copying
	 * @param offsetColumn column number to start copying
	 * @param row number of rows to be coppied
	 * @param column number of columns to be coppied
	 * @return the specified sub-matrix
	 */
	public IntMatrix getSubmatrix(int offsetRow, int offsetColumn, int row, int column){
		IntMatrix ret;
		
		ret = new IntMatrix(row, column);
		
		//copy the specified part of the matrix
		for(int i=0; i<row; i++){
			for(int j=0; j<column; j++){
				ret.data[i][j] = data[i+offsetRow][j+offsetColumn];
			}
		}
		
		return ret;
	}
	
	/**
	 * Calculates the maximum element
	 * @return the maximum element
	 */
	public int maxElement(){
		int maxEntry;
		maxEntry = this.data[0][0];
		for (int i = 0; i < this.row; i++){
			for (int j = 0; j < this.column; j++){
				if (maxEntry < this.data[i][j]){
					maxEntry = this.data[i][j];
				}
			}
		}
		return maxEntry;
	}
	
	/**
	 * Calculates the minimum element
	 * @return the minimum element
	 */
	public int minElement(){
		int minEntry;
		minEntry = this.data[0][0];
		for (int i = 0; i < this.row; i++){
			for (int j = 0; j < this.column; j++){
				if (minEntry > this.data[i][j]){
					minEntry = this.data[i][j];
				}
			}
		}
		
		return minEntry;
	}
	
	/**
	 * Produces a new matrix, transpose of the original one such that
	 * transposeMatrix[j][i] = originalMatrix[i][j] for all i,j 
	 * @return The transpose matrix
	 */
	public IntMatrix transpose(){
		IntMatrix trans = new IntMatrix(this.column, this.row);
		int[][] transData = trans.data;
		for(int i = 0; i<this.row; i++){
			for(int j=0; j<this.column; j++){
				transData[j][i] = data[i][j];
			}
		}
		
		return trans;
	}

	/**
	 * Replaces elements of the matrix which equals to formerData
	 * and sets them to newData
	 * @param formerData The element to be replaced
	 * @param newData The new element to be set
	 */
	public void setElement(int formerData, int newData){
		for(int i = 0; i<this.row; i++){
			for(int j=0; j<this.column; j++){
				if(data[i][j] == formerData){
					data[i][j] = newData;
				}
			}
		}
	}
	
	public IntMatrix submatrix(int startX, int startY){
		int crow = this.row - startX;
		int ccolumn = this.column - startY;
		
		IntMatrix ret = new IntMatrix(crow, ccolumn);
		
		for(int priX=startX, secX=0; priX<this.row; priX++, secX++)
			for(int priY=startY, secY=0; priY<this.column; priY++, secY++)
				ret.data[secX][secY] = this.data[priX][priY];
		
		return ret;
	}
	
	/**
	 * Reads a given file and loads the data into the data member.<br/>
	 * File format is as follows:<br/>
	 * Number_of_Rows Number_of_Columns<br/>
	 * x1 x2 x3.... (first row of data)
	 * @param filename File to be loaded.
	 * @return Itself. If an exception occurs (i.e. FileNotFoundException) 
	 * returns null.
	 */
	public static IntMatrix readFromFile(String filename){
		try{
			IntMatrix ret;
			FileReader fr = new FileReader(filename); 
			BufferedReader br = new BufferedReader(fr);
			String line = null;
			StringTokenizer tok = null;
			
			line=br.readLine();
			tok = new StringTokenizer(line);
			
			int crow = Integer.parseInt(tok.nextToken());
			int ccolumn = Integer.parseInt(tok.nextToken());
			ret = new IntMatrix(crow, ccolumn);
			
			for(int i=0; i < ret.row; i++){
				line = br.readLine();
				tok = new StringTokenizer(line);
				for(int j=0; j < ret.column; j++){
					ret.data[i][j] = Integer.parseInt(tok.nextToken());
				}
			}
			br.close();
			fr.close();
			return ret;
		}
		catch(java.io.FileNotFoundException e){
			e.printStackTrace();
			return null;
		}
		catch(java.io.IOException e){
			e.printStackTrace();
			return null;
		}
	}
	
	public static IntMatrix readFromHeaderlessFile(String filename, int crow, int ccolumn){
		try{
			IntMatrix ret;
			FileReader fr = new FileReader(filename); 
			BufferedReader br = new BufferedReader(fr);
			String line = null;
			StringTokenizer tok = null;
			
			ret = new IntMatrix(crow, ccolumn);
			
			for(int i=0; i < ret.row; i++){
				line = br.readLine();
				tok = new StringTokenizer(line);
				for(int j=0; j < ret.column; j++){
					ret.data[i][j] = Integer.parseInt(tok.nextToken());
				}
			}
			br.close();
			fr.close();
			return ret;
		}
		catch(java.io.FileNotFoundException e){
			e.printStackTrace();
			return null;
		}
		catch(java.io.IOException e){
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * Vertically concatenates two matrices.
	 * @param A First matrix
	 * @param B Second matrix
	 * @return vertically concatenated matrix
	 */
	public static IntMatrix vertcat(IntMatrix A, IntMatrix B){
		if(A.column != B.column){
			System.err.println("Matrix dimensions mismatch.");
			return null;
		}
		else{
			IntMatrix ret = new IntMatrix( (A.row + B.row), A.column);
			ret.copyMatrix(A, 0, 0);
			ret.copyMatrix(B, A.row, 0);
			return ret;
		}
	}
	
	/**
	 * Writes the data of the matrix to a specified file. <br/>
	 * File format is as follows: <br/>
	 * Number_of_Rows Number_of_Columns <br/>
	 * x1 x2 x3.... (first row of data)
	 * @param filename File name to be written
	 * @param matrixSize Whether the height and width of the matrix is going to be written or not
	 */
	public void writeToFile(String filename, boolean matrixSize){
		try{
			FileWriter fw = new FileWriter(filename); 
			BufferedWriter bw = new BufferedWriter(fw);
			StringBuffer sb = new StringBuffer();
			
			if(matrixSize){
				bw.write(row + " " + column + "\r\n");
			}
			
			for(int i = 0; i<this.row; i++){
				sb = new StringBuffer();
				sb.append(data[i][0]);
				//write a single column at once
				for(int j=1; j<this.column;j++){
					sb.append(" ");
					sb.append(data[i][j]);
				}
				sb.append("\r\n");
				bw.write(sb.toString());
			}
				
			bw.close();
			fw.close();
		}
		catch(java.io.IOException e){
			e.printStackTrace();
		}
	}

	public String toString(){
		String ret = "";
		ret += "Number of rows: " + row + "\n";
		ret += "Number of columns: " + column + "\n";
		
		for(int i = 0; i<this.row; i++){
			for(int j=0; j<this.column;j++){
				ret += data[i][j] + " ";
			}
			ret += "\n";
		}
		return ret;
	}
	
	/**
	 * Allocates the data member with given size 
	 * @param row
	 * @param column
	 */
	private void init(int row, int column){
		data = new int[row][column];
		this.row = row;
		this.column = column;
	}
}
