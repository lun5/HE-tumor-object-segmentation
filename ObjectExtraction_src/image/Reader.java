package image;

import java.awt.Image;
import java.awt.Toolkit;
import java.awt.image.PixelGrabber;


public class Reader {

	static HistoImage hImage=null;


	public static HistoImage readImage(String fileName, String d, int imageType){
		Image img = Toolkit.getDefaultToolkit().createImage(fileName);


		if (img!=null) {
			try {
				hImage = new HistoImage(fileName, img);
			} catch (IllegalStateException e) {
				
				return null; // error loading image
			}	
			catch(ImageDoesNotExistException e){
				System.out.println("Error loading image file.");
				System.exit(0);
			}
			int width = hImage.getWidth();
			int height = hImage.getHeight();
			int[] pixels =  new int[width * height];
			PixelGrabber pg = new PixelGrabber(img, 0, 0, width, height, pixels, 0, width);
			try {
				pg.grabPixels();
			} catch (InterruptedException e){};

			pixels = (int[])pg.getPixels();
			int imagesize = width * height;


			int x=0,y=0;
			for (int j=0;j<imagesize;j++){
				if(j>0 && j%width==0){ 
					y=0; x++;		
				}
				double R = (pixels[j] & 0xff0000)>>16;
				double G = (pixels[j] & 0x00ff00)>>8 ;
				double B = (pixels[j] & 0x0000ff);

				hImage.setR(x, y, R);
				hImage.setG(x, y, G);
				hImage.setB(x, y, B);
				y++;

			}
			hImage.setColorImage();

		}

		return hImage;
	}
}
