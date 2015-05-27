package main;



import image.HistoImage;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import util.Utility;

import matrix.Matrix;

public class Kmeans {



	final static int KMEANS_MAX_ITER_NO = 500;
	final static int RANDOM = 0;
	final static int MEAN = 1;
	final static int PCA = 2;
	
	public static void kmeansData(String dataFile, String mapFile, String vFile, int K, int initializationMethod){
		
		Matrix data=null;
		
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(dataFile));
			data				= Matrix.read(br,1);

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		Matrix V			= kmeans(data,K,initializationMethod);
		Matrix classMap		= findClusterLabelsUsingEuclideanDistance(data,V);
		
		classMap.writeMatrixIntoFile(mapFile,1,1);
		V.writeMatrixIntoFile(vFile,1,0);

		
	}
	
	public static Matrix kmeansImageOnTheFly(HistoImage hImage , String vFile, int K, int initializationMethod){
		Matrix gaborImage=null;
		Matrix data =null;
		int lm_id, st_id, nc_id, i, j;
		for(int k=1; k<=3; ++k){
				
				gaborImage	= hImage.getChannel(k);

			if(data==null){
				data = new Matrix(gaborImage.getRowDimension()*gaborImage.getColumnDimension(), 3, 0);
			}
			int count =0;
			for (i = 0; i < gaborImage.getRowDimension(); i++)
				for (j = 0; j < gaborImage.getColumnDimension(); j++){
					data.set(count, k-1, gaborImage.get(i, j));
					count++;
				}
		
		}
		
		Matrix V			= kmeans(data,K,initializationMethod);	
		Matrix classMap		= findClusterLabelsUsingEuclideanDistance(data,V);
		Matrix classMap2D	= classMap.form2DMatrix(gaborImage.getRowDimension(),gaborImage.getColumnDimension());

		
		if ((V.get(0,0) > V.get(1,0)) && (V.get(0,0) > V.get(2,0))){
			lm_id = 0;
			if (V.get(1,0) > V.get(2,0)){	st_id = 1;		nc_id = 2;	}
			else{								st_id = 2;		nc_id = 1;	}
		}
		else if ((V.get(1,0) > V.get(0,0)) && (V.get(1,0) > V.get(2,0))){
			lm_id = 1;
			if (V.get(0,0) > V.get(2,0)){	st_id = 0;		nc_id = 2;	}
			else{								st_id = 2;		nc_id = 0;	}
		}
		else{
			lm_id = 2;
			if (V.get(0,0) > V.get(1,0)){	st_id = 0;		nc_id = 1;	}
			else{								st_id = 1;		nc_id = 0;	}
		}


		for (i = 0; i < classMap2D.getRowDimension(); i++)
			for (j = 0; j < classMap2D.getColumnDimension(); j++){
				if (classMap2D.get(i, j) == lm_id)
					classMap2D.set(i, j, Utility.LM_PIX);
				else if (classMap2D.get(i, j) == nc_id)
					classMap2D.set(i, j, Utility.NC_PIX);
				else if (classMap2D.get(i, j) == st_id)
					classMap2D.set(i, j, Utility.ST_PIX);
				else{
					System.out.println("Error in getClusteredPixels (2)");
					System.exit(1);
				}
			}
		
		// YOU CAN COMMENT OUT THE FOLLOWING LINE TO NOT WRITING THE CLUSTERED PIXEL FILE
		classMap2D.writeMatrixIntoFile(vFile, 1, 1);
		return classMap2D;

	}
	
	public static void kmeansImage(String imFile, String mapFile, String vFile, int K, int initializationMethod){
		Matrix gaborImage=null;
		Matrix data =null;
		int i,j;
		for(int k=1; k<=12; ++k){
			BufferedReader br;
			try {
				br = new BufferedReader(new FileReader(imFile+"."+k));
				gaborImage				= Matrix.read(br,1);

			} catch (FileNotFoundException e) {
				e.printStackTrace();
				System.exit(1);
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}
			if(data==null){
				data = new Matrix(gaborImage.getRowDimension()*gaborImage.getColumnDimension(), 12, 0);
			}
			int count =0;
			for (i = 0; i < gaborImage.getRowDimension(); i++)
				for (j = 0; j < gaborImage.getColumnDimension(); j++){
					data.set(count, k-1, gaborImage.get(i, j));
					count++;
				}
		
		}
		
		Matrix V			= kmeans(data,K,initializationMethod);	
		Matrix classMap		= findClusterLabelsUsingEuclideanDistance(data,V);
		Matrix classMap2D	= classMap.form2DMatrix(gaborImage.getRowDimension(),gaborImage.getColumnDimension());

		
		classMap2D.writeMatrixIntoFile(mapFile, 1, 1);
		

	}
	public static Matrix kmeans(Matrix data, int K, int initializationMethod){
		Matrix B;
		Matrix M, prevM;
		double maxDeltaM;
		int iterNo;

		M = initializeClusterVectors(data,K,initializationMethod);
		prevM = new Matrix(K,data.getColumnDimension());
		
		
		iterNo = 0;
		while (iterNo <= KMEANS_MAX_ITER_NO){
			iterNo++;
			prevM = M.copy();
			B = new Matrix(data.getRowDimension(),K,0);
			estimateClustersUsingEuclideanDistance(B,data,M);
			computeClusterVectors(M,data,B);
			maxDeltaM = M.maxAbsoluteDifferenceBetweenMatrices(prevM);
			if (maxDeltaM <= Utility.SMALL)
				break;
		}
		
		return M;
	}
	public static void estimateClustersUsingEuclideanDistance(Matrix clusterAssigned, Matrix data, Matrix V){
		int i, j, t, minIndex, K = V.getRowDimension();
		double [] distance = new double[V.getRowDimension()];

		
		for (i = 0; i < data.getRowDimension(); i++){
			for (j = 0; j < K; j++){
				distance[j] = 0.0;
				for (t = 0; t < data.getColumnDimension(); t++)
					distance[j] += Utility.SQRDIST(V.get(j, t),data.get(i, t));
			}
			minIndex = Utility.indexOfMinArrayEntry(distance, K);
			clusterAssigned.set(i, minIndex, 1);
		}
	}
	public static Matrix findClusterLabelsUsingEuclideanDistance(Matrix data, Matrix V){
		int i, j, t, minIndex, K = V.getRowDimension();
		double [] distance = new double[V.getRowDimension()];
		Matrix labels = new Matrix(data.getRowDimension(),1);

		for (i = 0; i < data.getRowDimension(); i++){
			for (j = 0; j < K; j++){
				distance[j] = 0.0;
				for (t = 0; t < data.getColumnDimension(); t++)
					distance[j] += Utility.SQRDIST(V.get(j, t),data.get(i,t));
			}
			minIndex = Utility.indexOfMinArrayEntry(distance, K);
			labels.set(i, 0, minIndex);
		}
		return labels;
	}
	public static void computeClusterVectors(Matrix M, Matrix data, Matrix B){
		int i, j, t;
		Matrix N = new Matrix(M.getRowDimension(),1,0);

		M.initializeMatrix(0.0);
		

		for (t = 0; t < B.getRowDimension(); t++)
			for (i = 0; i < B.getColumnDimension(); i++)
				if ((int)B.get(t, i) == 1){
					for (j = 0; j < M.getColumnDimension(); j++)
						M.set(i, j, M.get(i, j)+data.get(t, j));//M->data[i][j] += data.data[t][j];
					N.set(i, 0, N.get(i, 0)+1);
				}
		for (i = 0; i < N.getRowDimension(); i++)
			if (N.get(i, 0) > 0)
				for (j = 0; j < M.getColumnDimension(); j++)
					M.set(i, j, M.get(i, j)/N.get(i, 0)); //M->data[i][j] /= N.data[i][0];
			else{		// randomly assign new clustering vector
				System.out.println("Randomly re-assigning clustering vectors...");
				
				t = (int)(Math.random() * data.getRowDimension());
				for (j = 0; j < M.getColumnDimension(); j++)
					M.set(i, j, data.get(t, j));//M->data[i][j] = data.data[t][j];
			}
	}
	public static Matrix initializeClusterVectors(Matrix data, int K, int initializationMethod){
		Matrix M=null;
		if (initializationMethod == RANDOM)
			M = initializeClusterVectorsRandomly(data,K);
		else if (initializationMethod == PCA)
			M = initializeClusterVectorsUsingPrincipalComponent(data,K);
		else{
			System.out.println("\nError: Invalid initialization method");
			System.exit(1);
		}
		return M;
	}
	public static Matrix initializeClusterVectorsRandomly(Matrix data, int K){
		int [] indices = new int [K];
		int i, j, t;
		Matrix M;
		
		
		for (i = 0; i < K; i++){
			while(true){
				int found = 0;
				t = (int)(Math.random() * data.getRowDimension());
				for (j = 0; j < i; j++)
					if (t == indices[j]){
						found = 1;
						break;
					}
				if (found == 0)
					break;
			}
			indices[i] = t;
		}
		// initialize cluster vectors V 
		M = new Matrix(K, data.getColumnDimension()); 
		for (i = 0; i < M.getRowDimension(); i++)
			for (j = 0; j < M.getColumnDimension(); j++)
				M.set(i, j, data.get(indices[i], j)); //M.data[i][j] = data.data[indices[i]][j];
		
		return M;
	}

	public static Matrix initializeClusterVectorsUsingPrincipalComponent(Matrix data, int K){
		Matrix M, corrData, V, newData, newM, temp;
		double d, maxEntry, minEntry, c[];
		int i, j, l, no[];
		
		corrData = data.computeCorrelationMatrix();
		V = new Matrix(corrData.getRowDimension(), corrData.getColumnDimension());

		corrData.computeEigenValues(V);
		newData = data.multiplyMatrix(V);
		maxEntry = newData.maxMatrixEntry(0);
		minEntry = newData.minMatrixEntry(0);
		d = (maxEntry - minEntry) / K;
		c = new double[K+1];
		c[0] = minEntry;
		for (i = 1; i <= K; i++)
			c[i] = c[i-1] + d;
		
		no = new int [K]; 
		newM = new Matrix (K,data.getColumnDimension(),0.0);
		
		for (i = 0; i < newData.getRowDimension(); i++)
			for (j = 0; j < K; j++)
				if ((newData.get(i,0) >= c[j]) && (newData.get(i,0) < c[j+1])){
					for (l = 0; l < newM.getColumnDimension(); l++)
						newM.set(j, l, newM.get(j,l) + newData.get(i,l));
					no[j]++;
					break;
				}

		for (i = 0; i < K; i++)
			for (j = 0; j < newM.getColumnDimension(); j++)
				newM.set(i, j, newM.get(i, j) / no[i]);
		temp = V.inverseMatrix();
		M = newM.multiplyMatrix(temp);
		

		return M;
	}

}
