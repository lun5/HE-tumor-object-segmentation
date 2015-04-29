package main;

import matrix.Matrix;
import util.Utility;

public class Circles {

	static Matrix map =null;		
	static Matrix currR =null;
	static Matrix marked =null;
	static int id =1;
	public static void locateOneTissueComponent(Matrix cMap,Matrix objects, Matrix objectMap,int componentLabel, 
			int radiusThr, int squareSize){

		marked	= new Matrix (cMap.getRowDimension(),cMap.getColumnDimension());
		Matrix temp, se;
		int i, j;


		for (i = 0; i < cMap.getRowDimension(); i++)
			for (j = 0; j < cMap.getColumnDimension(); j++)
				if (cMap.get(i, j) == componentLabel)// || cMap.get(i, j) == 1)
					marked.set(i, j, 1);
				else
					marked.set(i, j, 0);

		if (squareSize > 0){
			se = new Matrix(squareSize,squareSize,1);

			bwdilate2(marked,se);

		}
		for (i = 0; i < objectMap.getRowDimension(); i++)
			for (j = 0; j < objectMap.getColumnDimension(); j++)
				if (objectMap.get(i, j)!= Utility.MYBACKGROUND)
					marked.set(i, j, 0);

		//marked.writeMatrixIntoFile("marked_"+componentLabel, 0, 1);
		
		temp = locateCircles(objects,marked,radiusThr,componentLabel);

		for (i = 0; i < temp.getRowDimension(); i++)
			for (j = 0; j < temp.getColumnDimension(); j++)
				if (temp.get(i, j)>0){
					objectMap.set(i, j, componentLabel);
				}


	}

	private static Matrix locateCircles(Matrix objects, Matrix marked, int radiusThr, int componentLabel) {
		map		= new Matrix (marked.getRowDimension(),marked.getColumnDimension());
		currR	= new Matrix (marked.getRowDimension(),marked.getColumnDimension());
		int maxR		= findMaxPossibleRadius(marked,currR) + 1;

		CircularBoundary B[]	= createCircularBoundaries(maxR);
		int cnt, i, x, y;


		findExactRadius(currR,B,maxR);
		map.initializeMatrix(0);

		cnt = 1;
		for (i = maxR - 1; i >= radiusThr; i--){
			for (x = 0; x < currR.getRowDimension(); x++)
				for (y = 0; y < currR.getColumnDimension(); y++)
					if (currR.get(x, y) == i){
						locateOneCircle(objects, currR,map,x,y,B,i,cnt);
						cnt++;
						// write the circle information into the file
						// centX   centY   type   radius
						//fprintf(fid,"%d\t%d\t%d\t%d\n",x,y,componentLabel,i);
					}
		}

		return map;	
	}

	private static void locateOneCircle(Matrix objects, Matrix radii, Matrix map2, int x, int y, CircularBoundary[] B, int R, int label) {
		int i, j, cx, cy, mx, my, r;
		int minx, miny, maxx, maxy;

		for (i = 0; i <= R; i++){
			for (j = 0; j < B[i].N; j++){
				cx = B[i].x[j];
				cy = B[i].y[j];

				map.set(x+cx, y+cy, label);
				map.set(x+cx, y-cy, label);
				map.set(x-cx, y+cy, label);
				map.set(x-cx, y-cy, label);

				objects.set(x+cx, y+cy, id);
				objects.set(x+cx, y-cy, id);
				objects.set(x-cx, y+cy, id);
				objects.set(x-cx, y-cy, id);
				
				currR.set(x+cx, y+cy, -1);
				currR.set(x+cx, y-cy, -1);
				currR.set(x-cx, y+cy, -1);
				currR.set(x-cx, y-cy, -1);
				
			}
			
		}
		id = id +1;
		minx = x - (2*R);		if (minx < 0)				minx = 0;
		maxx = x + (2*R);		if (maxx >= map.getRowDimension())		maxx = map.getRowDimension() - 1;
		miny = y - (2*R);		if (miny < 0)				miny = 0;
		maxy = y + (2*R);		if (maxy >= map.getColumnDimension())	maxy = map.getColumnDimension() - 1;
		
		for (mx = minx; mx <= maxx; mx++)
			for (my = miny; my <= maxy; my++){
				if (currR.get(mx, my) < 0)
					continue;
				r = (int)currR.get(mx, my);
				for (i = 0; i <= r; i++){
					for (j = 0; j < B[i].N; j++){
						cx = B[i].x[j];
						cy = B[i].y[j];

						if ((currR.get(mx+cx, my+cy)== -1) || (currR.get(mx+cx, my-cy) == -1) ||
								(currR.get(mx-cx, my+cy) == -1) || (currR.get(mx-cx, my-cy) == -1))
							break;
					}
					if (j < B[i].N)
						break;
				}
				currR.set(mx, my, i-1);
			}
		
	}

	private static void findExactRadius(Matrix out, CircularBoundary[] B, int maxR) {
		int x, y, k, i, r, cx, cy;

		for (x = 0; x < currR.getRowDimension(); x++)
			for (y = 0; y < currR.getColumnDimension(); y++)
				if (currR.get(x, y) > 1){
					r = (int)currR.get(x, y);
					for (k = 2; k <= r; k++){
						for (i = 0; i < B[k].N; i++){
							cx = B[k].x[i];
							cy = B[k].y[i];
							if ((currR.get(x+cx, y+cy)== -1) || (currR.get(x+cx, y-cy) == -1) ||
								(currR.get(x-cx, y+cy) == -1) || (currR.get(x-cx, y-cy) == -1))
								break;
						}
						if (i < B[k].N)
							break;
					}
					currR.set(x, y, k-1);
				}
		
	}

	private static CircularBoundary[] createCircularBoundaries(int maxR) {

		CircularBoundary[] B = new CircularBoundary [maxR];
		Matrix  M = new Matrix(maxR,maxR,-1);
		int i, j, k, cnt;



		M.set(0, 0, 0);
		B[0]= new CircularBoundary();
		B[0].N = 1;
		B[0].r = 0;
		B[0].x = new int[1];
		B[0].y = new int[1];
		B[0].x[0] = 0;
		B[0].y[0] = 0;

		for (k = 1; k < maxR; k++){
			B[k]= new CircularBoundary();
			B[k].r = k;
			B[k].N = 0;
			B[k].x = B[k].y = null;
			for (i = 0; i <= k; i++)
				for (j = 0; j <= k; j++){
					if (M.get(i, j) != -1)
						continue;
					if (i*i + j*j <= k*k){
						M.set(i, j, k);
						B[k].N += 1;
					}
				}
			if (B[k].N == 0)
				continue;
			B[k].x = new int[B[k].N];
			B[k].y = new int[B[k].N];
			cnt = 0;
			for (i = 0; i <= k; i++)
				for (j = 0; j <= k; j++)
					if (M.get(i, j)== k){
						B[k].x[cnt] = i;
						B[k].y[cnt] = j;
						cnt++;
					}
		}
		//writeMatrixIntoFile(M,"temp",0);

		return B;
	}


	private static int findMaxPossibleRadius(Matrix marked, Matrix out) {
		int maxR = 0, i, j, k;

		currR.initializeMatrix(-1);
		for (i = 0; i < marked.getRowDimension(); i++)
			for (j = 0; j < marked.getColumnDimension(); j++){
				if (marked.get(i, j) == 0)
					continue;

				k = 1;
				while(true){
					if ((i-k < 0)	|| (marked.get(i-k, j)== 0))		
						break;
					if ((i+k >= marked.getRowDimension())	|| (marked.get(i+k,j) == 0))		
						break;
					if ((j-k < 0)	|| (marked.get(i,j-k) == 0))		
						break;
					if ((j+k >= marked.getColumnDimension())|| (marked.get(i,j+k) == 0))		
						break;
					k++;
				}
				if(k>50) //you can set upper limit for the radius of the max possible radius
					k=50;
				currR.set(i, j, k-1);
				if (maxR < k-1)
					maxR = k-1;
			}
		return maxR;
	}

	private static void bwdilate2(Matrix marked2, Matrix se) {
		Matrix R = new Matrix(marked.getRowDimension(),marked.getColumnDimension(),0);
		int i, j, k1, k2;
		int msize = se.getRowDimension() / 2;

		if (se.getRowDimension() % 2 == 0){	// if even, exit
			System.out.println("\nError: Size of the structural element should be odd\n\n");
			System.exit(1);
		}

		for (i = 0; i < marked.getRowDimension(); i++)
			for (j = 0; j < marked.getColumnDimension(); j++){
				int found = 0;
				if (marked.get(i, j) != 0 && marked.get(i, j) != 1){
					System.out.println("\nError: Only BW images can be dilated\n\n");
					System.exit(1);
				}
				for (k1 = -msize; k1 <= msize; k1++){
					for (k2 = -msize; k2 <= msize; k2++){
						int x = i + k1; 
						int y = j + k2;
						if ((x < 0) || (x >= marked.getRowDimension()) || (y < 0) || (y >= marked.getColumnDimension()))
							continue;
						if (se.get(k1 + msize, k2 + msize)>0 && (marked.get(x, y) == 1)){
							R.set(i, j, 1);
							found = 1;
							break;
						}
					}
					if (found==1)
						break;
				}
			}
		for (i = 0; i < marked.getRowDimension(); i++)
			for (j = 0; j < marked.getColumnDimension(); j++)
				marked.set(i, j, R.get(i, j));


	}



}
