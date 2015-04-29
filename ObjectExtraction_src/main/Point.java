package main;

public class Point {

	int X, Y, rNo;
	int visited=0;
	public Point(int x, int y, int no) {
		super();
		X = x;
		Y = y;
		rNo = no;
	}

	public int getRNo() {
		return rNo;
	}

	public void setRNo(int no) {
		rNo = no;
	}

	public int getX() {
		return X;
	}

	public void setX(int x) {
		X = x;
	}

	public int getY() {
		return Y;
	}

	public void setY(int y) {
		Y = y;
	}
	
	
}
