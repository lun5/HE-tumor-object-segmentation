package image;



import java.awt.Canvas;
import java.awt.Color;
import java.awt.Component;
import java.awt.Frame;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;

import matrix.Matrix;

public class HistoImage implements ImageObserver{

	protected Image img;
	public Image getImg() {
		return img;
	}

	public void setImg(Image img) {
		this.img = img;
	}

	private String title;
	protected int width;
	protected int height;
	private int ID;
	private Matrix R=null;
	private Matrix G=null;
	private Matrix B=null;
	private int colorImage[][];
	private static Component comp;
	
	private static int currentID = -1;
	
	public HistoImage(String title, Image img) throws ImageDoesNotExistException {
		this.title = title;
    	ID = --currentID;
    	int check=1;
		if (img!=null)
			check=setImage(img);
		if(check==-1)
			throw new ImageDoesNotExistException();
		colorImage= new int[height][width];
    }
	
	public int setImage(Image img) {
		comp = new Canvas();
		int count =0;
		while (!comp.prepareImage(img, this)) {
			count++;
			if(count> 900000000)
				return -1;
		}
		this.img = img;
		int newWidth = img.getWidth(this);
		int newHeight = img.getHeight(this);
		boolean dimensionsChanged = newWidth!=width || newHeight!=height;
		width = newWidth;
		height = newHeight;
		R= new Matrix(height,width);
		G= new Matrix(height,width);
		B= new Matrix(height,width);
		return 1;
	}

	public int getHeight() {
		return height;
	}

	public void setHeight(int height) {
		this.height = height;
	}

	public int getWidth() {
		return width;
	}

	public void setWidth(int width) {
		this.width = width;
	}
	 public boolean imageUpdate(Image img, int flags, int x, int y, int w, int h) {

			
			return true;
	    }
	 public void setR (int i, int j , double a){
		 R.set(i, j, a);
	 }
	 public void setG (int i, int j , double a){
		 G.set(i, j, a);
	 }
	 public void setB (int i, int j , double a){
		 B.set(i, j, a);
	 }

	 public Matrix getChannel(int channel){
		 if(channel==1)
			 return R;
		 else if(channel==2)
			 return G;
		 else
			 return B;
		 
	 }
	public Matrix getB() {
		return B;
	}

	public void setB(Matrix b) {
		B = b;
	}

	public Matrix getG() {
		return G;
	}

	public void setG(Matrix g) {
		G = g;
	}

	public Matrix getR() {
		return R;
	}

	public void setR(Matrix r) {
		R = r;
	}

	public int[][] getColorImage() {
		return colorImage;
	}

	public void setColorImage() {
		for(int i=0; i< R.getRowDimension(); i++){
			for(int j=0; j< R.getColumnDimension(); j++){
				Color c= new Color((int)R.get(i,j) , (int)G.get(i,j), (int)B.get(i,j));
				colorImage[i][j]=c.getRGB();
			}
		}
	}
	

}
