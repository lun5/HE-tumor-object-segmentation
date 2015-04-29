package main;

import image.HistoImage;
import image.Reader;

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Vector;

import javax.imageio.ImageIO;






import matrix.*;
import feature.*;
import util.*;
import delaunay.*;

public class Runlength {

	/**
	 * @param args
	 * @return 
	 */


	static Matrix cmap=null;
	static Matrix objects=null;
	static Matrix goldStandart=null;
	static ImageObject objProp[]=null;
	static ImageObject superObjProp[]=null;

	static Matrix objectsCreated=null;
	static Matrix objectMap=null;

	public static void main(String[] args) {



		if (args.length<6){
			System.out.println("Need more parameters.");
			System.out.println("Parameters: <jpeg filename without extension> <min. nuclei radius> <min. stroma radius> <min. lumen radius> <struct element radius>  <window size> ");
			System.exit(0);
		}	



		int k=3;
		String fileName = args[0]; 
		String imageFile = fileName +".jpg";
		int minNucleiRadius = Integer.parseInt(args[1]);
		int minStromaRadius = Integer.parseInt(args[2]);
		int minLumenRadius = Integer.parseInt(args[3]);
		int structElementRadius = Integer.parseInt(args[4]);
		String kmapFile = fileName + "_k" + Integer.toString(k);
		fileName = fileName+"_se"+structElementRadius+"_minNuc"+minNucleiRadius+"_minStr"+minStromaRadius+"_minLum"+minLumenRadius;
		Calendar c1 = Calendar.getInstance();
		//
		String cmapFile = fileName+ "_circle_map";
		String omapFile = fileName + "_objects";

		Matrix cMap= null;
		File pf= new File(cmapFile);
		HistoImage h = Reader.readImage(imageFile, "rgb", Utility.HE);

		if(!pf.exists()){

			File pf2 = new File(kmapFile);

			if(!pf2.exists()){
				cMap = Kmeans.kmeansImageOnTheFly(h, kmapFile, k, Kmeans.PCA);
			}
			else{
				try {
					BufferedReader br = new BufferedReader(new FileReader(kmapFile));
					cMap=  Matrix.read(br,1);
				} catch (FileNotFoundException e) {
					e.printStackTrace();
					System.exit(1);
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(1);
				}
			}
			System.out.println("Starting object extraction...");
			createCircularObjects(fileName,cMap,minNucleiRadius,minStromaRadius,minLumenRadius,structElementRadius);
		}

		Matrix voronoiMap =null;
		Calendar c6 = Calendar.getInstance(); 
		long objectExtractTime = c6.getTimeInMillis()-c1.getTimeInMillis();
		System.out.println("Object exraction for image : " + fileName + " completed in " + objectExtractTime/1000 + " seconds");
		//System.exit(1);
		readFilesCmap(fileName);

		System.out.println("Computing object properties");
		voronoiMap=createObjectsWithPropsCmap(fileName, voronoiMap);
		System.out.println("Writing voronoi map  of objects into file.");
		voronoiMap.writeMatrixIntoFile(fileName+"_voronoiMap", 0,1);

		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(cmapFile));
			cmap				= Matrix.read(br,1);
			br = new BufferedReader(new FileReader(omapFile));
			objects				= Matrix.read(br,1);
  
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}

		objects.relabelComponents();

		Calendar c2 = Calendar.getInstance(); 
		long readWriteTime = c2.getTimeInMillis()-c1.getTimeInMillis();
		int onlyFeatureExtraction= 0 ; //Integer.parseInt(args[6]);

		int winSize = Integer.parseInt(args[5]);
		double distThr  = 1.25;
		int compThr = 100;
		System.out.println("Computing object graph");
		Matrix res =CIRCLE_REGION_GROWING_WITH_OBJECT_RUNLENGTH_FEATURES(fileName,cmap,objects, voronoiMap,
				distThr , winSize,compThr,onlyFeatureExtraction);

		//remove comments below to write result label map
		//res.writeMatrixIntoFile(fileName+"_result", 0,1);


		/*int noOfRegions = res.maxMatrixEntry();
		int totalRegionNorm[] = new int[noOfRegions+1];
		int totalRegionCanc[] = new int[noOfRegions+1];
		for (int ii = 0; ii < res.getRowDimension(); ii++)
			for (int jj = 0; jj < res.getColumnDimension(); jj++){
				int rno= (int)res.get(ii, jj);
				if(goldStandart.get(ii, jj)==1)
					totalRegionCanc[rno]++;
				if(goldStandart.get(ii, jj)==2)
					totalRegionNorm[rno]++;
			}

		int actualRegionType[] = new int[noOfRegions+1];
		for (int ii = 1; ii <= noOfRegions; ii++){
			if(totalRegionCanc[ii]>=totalRegionNorm[ii])
				actualRegionType[ii]=1;
			else
				actualRegionType[ii]=2;
		}

		double TP=0, TN=0, FP=0, FN=0;
		for (int ii = 0; ii < res.getRowDimension(); ii++)
			for (int jj = 0; jj < res.getColumnDimension(); jj++){
				int rno= (int)res.get(ii, jj);
				if(goldStandart.get(ii, jj)==1 && actualRegionType[rno]==1)
					TP++;
				if(goldStandart.get(ii, jj)==1 && actualRegionType[rno]==2)
					FN++;
				if(goldStandart.get(ii, jj)==2 && actualRegionType[rno]==1)
					FP++;
				if(goldStandart.get(ii, jj)==2 && actualRegionType[rno]==2)
					TN++;
			}

		double accuracy=0, sensitivity=0, specificity=0, dice =0, PPV=0 , PLR =0;

		accuracy= (TP+TN)/(TP+TN+FP+FN);
		sensitivity= TP/(TP+FN);
		specificity= TN/(TN+FP);

		DecimalFormat df = new DecimalFormat("#.####"); */

		//UNCOMMENT FOLLOWING COMMENTED LINES IF YOU HAVE GOLDSTANDART TO COMPUTE QUANTITATIVE RESULTS

		/*PrintWriter out=null;
		try {


			out = new PrintWriter(new FileWriter(fileName+"_batchResults_"+args[1],true));
			out.print("id="+fileName+",distThr="+Double.parseDouble(args[3])+",areaThr="+Integer.parseInt(args[5])+",winSize="+Integer.parseInt(args[4])+
				",growRate="+0.1+",");		
			out.print("method="+1+",TP="+TP+",FP="+FP+",");
			out.print("TN="+TN+",FN="+FN+"\t");

			out.print(df.format(accuracy)+"\t");

			if ((TP==0 &&  FN==0))
				out.print("0\t");
			else
				out.print(df.format(sensitivity)+"\t");
			if ((TN==0 &&  FP==0))
				out.print("0\t");
			else
				out.print(df.format(specificity)+"\t");

			if ((TP==0 &&  FN==0 &&  FP==0))
				out.print("0\t");
			else
				out.print(df.format(dice)+"\t");

			if(TP==0 && FP==0)
				out.print("0\t");
			else
				out.print(df.format(PPV)+"\t");	
			if(TN==0 && FN==0)
				out.println("1\t"+noOfRegions + "\t"+1 );
			else
				out.println(df.format(PLR)+ "\t"+noOfRegions + "\t"+1 );

			out.close();
		} catch (IOException e) {

			e.printStackTrace();
		}*/

	/*	BufferedImage theImage = new BufferedImage(objects.getColumnDimension(), objects.getRowDimension(), BufferedImage.TYPE_INT_RGB);
		for(int y = 0; y<objects.getRowDimension(); y++){
			for(int x = 0; x<objects.getColumnDimension(); x++){
				int value = (int)(h.getR().get(y,x)*h.getG().get(y,x)*h.getB().get(y,x));
				theImage.setRGB(x, y, h.getColorImage()[y][x]);
			}
		}

		Matrix edges = res.findEdges();
		edges.bwdilate2(new Matrix(9,9,1));
		for(int y = 0; y<objects.getRowDimension(); y++){
			for(int x = 0; x<objects.getColumnDimension(); x++){
				if(edges.get(y,x)==1)
					theImage.setRGB(x, y, 0);
			}
		}
		File outputfile = new File(fileName+"_win"+ winSize + "_distThr"+ distThr + "_compThr"+ compThr+"_segmentation_result.jpg");
		try {
			ImageIO.write(theImage, "jpg", outputfile);
		} catch (IOException e) {
			e.printStackTrace();
		}
*/
		Calendar c3 = Calendar.getInstance();
		long time = c3.getTimeInMillis()-c1.getTimeInMillis();
		System.out.println("Overall execution (with file read/write) took: "+ (time+readWriteTime)/1000+" secs.");
	}

	public static void createCircularObjects(String fileName, Matrix cMap, int minNucleiRadius, int minStromaRadius, int minLumenRadius, int structElementRadius) {
		objectsCreated	= new Matrix(cMap.getRowDimension(),cMap.getColumnDimension(),0);
		objectMap	= new Matrix(cMap.getRowDimension(),cMap.getColumnDimension(),Utility.MYBACKGROUND);

		System.out.println("Extracting nuclei objects.");
		Circles.locateOneTissueComponent(cMap,objectsCreated,objectMap,Utility.NC_PIX,minNucleiRadius,structElementRadius*2+1);
		System.out.println("Extracting nuclei objects.Done.");
		System.out.println("Extracting stroma objects.");
		Circles.locateOneTissueComponent(cMap,objectsCreated,objectMap,Utility.ST_PIX,minStromaRadius,structElementRadius*2+1);
		System.out.println("Extracting stroma objects.Done");
		System.out.println("Extracting lumen objects.");
		Circles.locateOneTissueComponent(cMap,objectsCreated,objectMap,Utility.LM_PIX,minLumenRadius,structElementRadius*2+1);
		System.out.println("Extracting lumen objects.Done.");
		System.out.println("Writing object maps into files...");
		objectsCreated.writeMatrixIntoFile(fileName+"_objects", 1, 1);
		objectMap.writeMatrixIntoFile(fileName+"_circle_map", 1, 1);

	}

	public static Matrix createObjectsWithPropsCmap(String fileName, Matrix vMap) {

		int i, j;

		int numOfNodes;

		objects.relabelComponents();
		objProp			= ImageObject.computeObjectPropertiesL(cmap,objects,objects.maxMatrixEntry());

		numOfNodes = objects.maxMatrixEntry();



		Matrix voronoiMap = new Matrix(objects.getRowDimension(), objects.getColumnDimension(),0);
		int arraySize = objects.getRowDimension()* objects.getColumnDimension();
		Point [] pointTail2 = new Point[arraySize];

		int countInit =0;
		for(i=1 ; i< numOfNodes+1 ;++i){
			ImageObject obj = objProp[i];
			obj.createNeighs(0);
			obj.noNeighbours=0;

			voronoiMap.set((int)obj.getCentroidX(), (int)obj.getCentroidY(), obj.regionNo);

			pointTail2[countInit] = new Point((int)obj.getCentroidX(), (int)obj.getCentroidY(), obj.regionNo);
			countInit++;
		}



		int countPoint = 0;
		if(vMap==null){

			while(countPoint < arraySize) { 

				Point p = pointTail2[countPoint];

				countPoint++;
				//four connectivity
				if(p.getX()+1 < voronoiMap.getRowDimension() ){
					if(voronoiMap.get(p.getX()+1, p.getY())==0){
						voronoiMap.set(p.getX()+1, p.getY(), p.getRNo());

						pointTail2[countInit] =new Point(p.getX()+1, p.getY(), p.getRNo());
						countInit++;
					}

				}
				if(p.getX()-1 >= 0  ){
					if(voronoiMap.get(p.getX()-1, p.getY())==0){
						voronoiMap.set(p.getX()-1, p.getY(), p.getRNo());

						pointTail2[countInit] =new Point(p.getX()-1, p.getY(), p.getRNo());
						countInit++;
					}
				}
				if(p.getY()+1 < voronoiMap.getColumnDimension() ){
					if(voronoiMap.get(p.getX(), p.getY()+1)==0){
						voronoiMap.set(p.getX(), p.getY()+1, p.getRNo());

						pointTail2[countInit] =new Point(p.getX(), p.getY()+1, p.getRNo());
						countInit++;
					}
				}
				if(p.getY()-1 >= 0  ){
					if(voronoiMap.get(p.getX(), p.getY()-1)==0){
						voronoiMap.set(p.getX(), p.getY()-1, p.getRNo());

						pointTail2[countInit] =new Point(p.getX(), p.getY()-1, p.getRNo());
						countInit++;
					}

				}
			}

		}
		if(vMap!=null){
			voronoiMap=vMap;
		}



		for (i = 0; i < objects.getRowDimension(); i++){
			objProp[(int)voronoiMap.get(i, 0)].isBorder=true;
		}
		for (i = 0; i < objects.getRowDimension(); i++){
			objProp[(int)voronoiMap.get(i, objects.getColumnDimension()-1)].isBorder=true;
		}
		for (j = 0; j < objects.getColumnDimension(); j++){
			objProp[(int)voronoiMap.get(0, j)].isBorder=true;
		}
		for (j = 0; j < objects.getColumnDimension(); j++){
			objProp[(int)voronoiMap.get(objects.getRowDimension()-1, j)].isBorder=true;
		}


		Triangle tri =
				new Triangle(new Pnt(-1,-1000000,1000000), new Pnt(-1,1000000,1000000), new Pnt(-1,0,-1000000));
		Triangulation dt = new Triangulation(tri);
		for(i=1 ; i< numOfNodes+1 ;++i){
			ImageObject obj = objProp[i];
			try{
				dt.delaunayPlace(new Pnt(obj.getRegionNo(),obj.getCentroidX(),obj.getCentroidY()));	
			}
			catch(Exception n) {

			}
		}
		i=1;

		for (Triangle triangle : dt) {
			Pnt[] polygon = triangle.toArray(new Pnt[0]);

			int cnt=0;
			if(polygon[1].getRNo()!=-1 && polygon[0].getRNo()!=-1 && polygon[2].getRNo()!=-1){
				if(!objProp[polygon[0].getRNo()].isBorder || !objProp[polygon[1].getRNo()].isBorder){
					objProp[polygon[0].getRNo()].updateNeighs(objProp[polygon[1].getRNo()]);
					objProp[polygon[1].getRNo()].updateNeighs(objProp[polygon[0].getRNo()]);
					cnt++;
				}

				if(!objProp[polygon[0].getRNo()].isBorder || !objProp[polygon[2].getRNo()].isBorder){
					objProp[polygon[0].getRNo()].updateNeighs(objProp[polygon[2].getRNo()]);
					objProp[polygon[2].getRNo()].updateNeighs(objProp[polygon[0].getRNo()]);
					cnt++;
				}


				if(!objProp[polygon[1].getRNo()].isBorder || !objProp[polygon[2].getRNo()].isBorder){
					objProp[polygon[2].getRNo()].updateNeighs(objProp[polygon[1].getRNo()]);
					objProp[polygon[1].getRNo()].updateNeighs(objProp[polygon[2].getRNo()]);
					cnt++;
				}
				if (cnt==3){

					i++;
				}
			}

		}

		return voronoiMap;
	}
	public static void readFilesCmap(String fileName1) {
		String fileName=fileName1;
		String cmapFile = fileName+ "_circle_map";
		String omapFile = fileName + "_objects";


		//UNCOMMENT FOLLOWING COMMENTED LINES IF YOU HAVE GOLDSTANDART FILE

		//String gsFile = fileName + "_gs";


		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(cmapFile));
			cmap				= Matrix.read(br,1);
			br = new BufferedReader(new FileReader(omapFile));
			objects				= Matrix.read(br,1);
			//br = new BufferedReader(new FileReader(gsFile));
			//goldStandart				= Matrix.read(br,1);

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}

		objects.relabelComponents();


	}

	static Matrix CIRCLE_REGION_GROWING_WITH_OBJECT_RUNLENGTH_FEATURES(String fileName, Matrix cmap, Matrix objects,  Matrix vMap,
			double thr, int winSize, int areaThre, int onlyFeatureExtraction){


		int i, j;
		int featNo =16;
		Matrix[] runLengthMatrices = new Matrix[7];
		int numOfNodes;
		BufferedReader br;

		objects.relabelComponents();

		numOfNodes = objects.maxMatrixEntry();

		Matrix featureM = null;

		File pf2 = new File(fileName+"_win"+winSize+"_runlength_features");
		System.out.println("Computing object graph runglength matrices and features.");
		if(!pf2.exists()){
			runLengthMatrices = GraphProp.findMaxRunLengthMatrixForEachNode(objProp,winSize);

			for( i=1; i<runLengthMatrices.length; i++){
				runLengthMatrices[i] = Matrix.deleteTheColumnsWithAllZeros(runLengthMatrices[i]); 
			}


			for(i=1 ; i< numOfNodes+1 ;++i){
				ImageObject obj = objProp[i];
				double [] featureVector = getObjFeatureVector(runLengthMatrices, objProp,  obj.getCentroidX(),obj.getCentroidY(), numOfNodes,winSize,featNo);

				obj.ptex= new double[featureVector.length];

				for (int k = 0; k < obj.ptex.length; k++)
					obj.ptex[k]= featureVector[k];
			}


			double [][] featureMatrix = new double[numOfNodes+1][objProp[1].ptex.length];

			for(i=1 ; i< numOfNodes+1 ;++i){
				ImageObject obj = objProp[i];
				for (int k = 0; k < obj.ptex.length; k++)
					featureMatrix[i][k]=obj.ptex[k];
			}

			featureMatrix = normalize(featureMatrix);

			featureM = new Matrix(numOfNodes,featNo);
			for(i=1 ; i< numOfNodes+1 ;++i){
				ImageObject obj = objProp[i];
				for (int k = 0; k < featNo; k++)
					featureM.set(i-1, k, featureMatrix[i][k]);
			}
			//featureM.writeMatrixIntoFile(fileName+"_win"+winSize+"_runlength_features", 1,0);
			if(onlyFeatureExtraction==1){		
				System.exit(0);
			}



			for(i=1 ; i< numOfNodes+1 ;++i){
				ImageObject obj = objProp[i];
				for (int k = 0; k < obj.ptex.length; k++)
					obj.ptex[k]= featureMatrix[i][k];
			}
		}
		else{
			try {
				br = new BufferedReader(new FileReader(fileName+"_win"+winSize+"_runlength_features"));
				featureM=  Matrix.read(br,1);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
				System.exit(1);
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}

			for(i=1 ; i< numOfNodes+1 ;++i){
				ImageObject obj = objProp[i];
				obj.ptex= new double[featureM.getColumnDimension()];
				for (int k = 0; k < featureM.getColumnDimension(); k++){

					obj.ptex[k]= featureM.get(i-1, k);
				}

			}
		}
		System.out.println("Computing object graph runglength matrices and features.Done.");
		System.out.println("Finding 50 closest neighbors and extracting their information..");
		Matrix neighX = new Matrix(objProp.length, 51, 0);
		Matrix neighY = new Matrix(objProp.length, 51, 0);
		Matrix neighFeatRun = new Matrix(objProp.length, 51, 0);
		Matrix neighRad = new Matrix(objProp.length, 51, 0);
		Matrix neighType = new Matrix(objProp.length, 51, 0);
		Matrix neighID = new Matrix(objProp.length, 51, 0);
		//for(i=1 ; i< numOfNodes+1 ;++i){
			//ImageObject obj = objProp[i];
			//neighs
			double[] distance;
			int[] nodes;
			int nodes_start = 1, nodes_end = 1, cur = 1;
			int countNodes=0;
			//for each connection type...
			//for(int connType=1; connType<=Utility.NUM_OF_CONN_TYPES; connType++){
			Matrix isUsed = new Matrix(objProp.length, 1, -1);

			//start from the first node
			for(int curRoot=1; curRoot<objProp.length; ++curRoot){
				countNodes=0;
				nodes = new int[objProp.length];
				distance = new double[52];
				//Vector<Integer> nodesInWindow = new  Vector<Integer>();
				//add the root of the tree (current node i) 
				//with a distance of zero
				nodes[1] = curRoot; distance[1] = 0;
				isUsed = new Matrix(objProp.length, 1, -1);
				isUsed.set(curRoot, 0, 1);
				nodes_start = 1; nodes_end = 1;
				double mainCenterX = objProp[curRoot].getCentroidX();
				double mainCenterY = objProp[curRoot].getCentroidY();

				for(;nodes_start<=nodes_end; nodes_start++){

					cur = nodes[nodes_start];
					boolean hasUnvisitedNeighbour = false;

					for(int k=0; k<objProp[cur].getNeigh().length; k++){

						double centX = objProp[objProp[cur].getNeigh()[k]].getCentroidX();
						double centY = objProp[objProp[cur].getNeigh()[k]].getCentroidY();
						double dist = Math.sqrt(Math.pow(centX-mainCenterX, 2) + Math.pow(centY-mainCenterY, 2));
						if(		(int)isUsed.get(objProp[cur].getNeigh()[k], 0) == -1 &&
								countNodes<50){
							//TODO inside window runlength control
							isUsed.set(objProp[cur].getNeigh()[k], 0, 0);

							//nodesInWindow.add(objProp[cur].getNeigh()[k]);
							hasUnvisitedNeighbour = true; 
							countNodes++;
							nodes[++nodes_end] = objProp[cur].getNeigh()[k];
							distance[nodes_end] = dist;

						}
					}

				}
				ArrayIndexComparator comparator = new ArrayIndexComparator(distance);
				Integer[] indexes = comparator.createIndexArray();
				Arrays.sort(indexes, comparator);
				if (countNodes!=0){
					countNodes=0;

					for (int jj=1;jj<indexes.length;jj++){
						neighX.set(curRoot, jj-1, objProp[nodes[indexes[jj]]].getCentroidX());
						neighY.set(curRoot, jj-1, objProp[nodes[indexes[jj]]].getCentroidY());
						double tempAre = objProp[nodes[indexes[jj]]].getArea();
						long tempR = Math.round(Math.sqrt(tempAre/Utility.PI));
						neighRad.set(curRoot, jj-1, tempR);
						neighFeatRun.set(curRoot, jj-1, objProp[nodes[indexes[jj]]].getPtex()[14]);
						neighType.set(curRoot, jj-1, objProp[nodes[indexes[jj]]].getClusterNo());
						neighID.set(curRoot, jj-1, objProp[nodes[indexes[jj]]].getObjectID());
					}
				}

			 

		}

		neighFeatRun.writeMatrixIntoFile(fileName+"_first50neigh_runlength_features",0,0);
		neighType.writeMatrixIntoFile(fileName+"_first50neigh_cluster_labels",0,1);
		neighRad.writeMatrixIntoFile(fileName+"_first50neigh_radii",0,0);
		neighX.writeMatrixIntoFile(fileName+"_first50neigh_Xcoors",0,1);
		neighY.writeMatrixIntoFile(fileName+"_first50neigh_Ycoors",0,1);
		neighID.writeMatrixIntoFile(fileName+"_first50neigh_IDs",0,1);
		System.out.println("Finding 50 closest neighbors and extrating their information. Done.");
		System.out.println("All Done for image : " + fileName );
		return  neighFeatRun;
		/*findSimlarityDistance(objProp, objects);


		/////////////


		Matrix seedMatrix = regionGrowing(objProp,objects,thr);



		seedMatrix.removeSmallNumberOfObjects(areaThre,objects);

		//You can see seed objects uncommenting next line
		//seedMatrix.writeMatrixIntoFile(fileName+"_seeds", 0,1);

		seedMatrix.relabelComponents();


		for (i = 0; i < objects.getRowDimension(); i++)
			for (j = 0; j < objects.getColumnDimension(); j++){
				int rno= (int)objects.get(i, j);
				int seedNo = (int)seedMatrix.get(i, j);
				if(rno>0 ){
					if(seedNo >0)
						objProp[rno].setNewRegionNo(seedNo);
					else
						objProp[rno].setNewRegionNo(-1);
				}
			}



		int maxSeedNo = seedMatrix.maxMatrixEntry();


		if (maxSeedNo <1)
			return new Matrix(objects.getRowDimension(),objects.getColumnDimension(),1.0);
		double [][] seedProp = new double[maxSeedNo+1][objProp[1].ptex.length];
		double [][] seedFirstProp = new double[maxSeedNo+1][objProp[1].ptex.length];
		int [] countSeedObjects = new int [maxSeedNo+1];
		int countSeedObj =0;
		for(i=1 ; i< numOfNodes+1 ;++i){
			ImageObject obj = objProp[i];
			if(obj.getNewRegionNo()>0){
				for (int k = 0; k < obj.ptex.length; k++){
					seedProp[obj.getNewRegionNo()][k]+= featureM.get(i-1, k);
				}
				countSeedObjects[obj.getNewRegionNo()]++;
				countSeedObj++;
			}
		}


		for(i=1; i<maxSeedNo+1; ++i){
			for (int k = 0; k < objProp[1].ptex.length; k++){
				seedFirstProp[i][k] = seedProp[i][k] / (double)countSeedObjects[i];
				seedProp[i][k] = seedProp[i][k] / (double)countSeedObjects[i];	
			}
		}


		double growThr = thr;
		double offset = growThr*0.5;
		int preNo =countSeedObj;
		while(countSeedObj<numOfNodes-1 && growThr<999){
			for(i=1 ; i< numOfNodes+1 ;++i){
				ImageObject obj = objProp[i];

				if(obj.getNewRegionNo()>0){			
					for(int k=0 ; k<obj.noNeighbours;++k){
						ImageObject obj2= objProp[obj.getNeigh()[k]];
						if(obj2.getNewRegionNo()<0){
							double d = 0, s=0;


							for (int x = 0; x < obj.ptex.length; x++)
								d += Utility.SQUARE(seedProp[obj.getNewRegionNo()][x] - obj2.ptex[x]);

							s = Math.sqrt(d);


							if(s < growThr){
								obj2.setNewRegionNo(obj.getNewRegionNo());
								for (int l = 0; l < objProp[1].ptex.length; l++){
									seedProp[obj.getNewRegionNo()][l] = seedProp[obj.getNewRegionNo()][l] * (double)countSeedObjects[obj.getNewRegionNo()];
									seedProp[obj.getNewRegionNo()][l] +=obj2.ptex[l];
									seedProp[obj.getNewRegionNo()][l] = seedProp[obj.getNewRegionNo()][l] / ((double)countSeedObjects[obj.getNewRegionNo()]+1.0);

								}
								countSeedObjects[obj.getNewRegionNo()]++;
								countSeedObj++;
							}

						}
					}
				}

			}

			if(preNo ==countSeedObj){
				growThr=growThr+offset;
			}
			else
				preNo =countSeedObj;
		}


		Matrix regMatrix = new Matrix(objects.getRowDimension(),objects.getColumnDimension(),0.0);
		Matrix vorMatrix = new Matrix(objects.getRowDimension(),objects.getColumnDimension(),0.0);
		for (i = 0; i < objects.getRowDimension(); i++)
			for (j = 0; j < objects.getColumnDimension(); j++){
				int rno= (int)objects.get(i, j);			
				if(rno>0 ){
					int seedNo = objProp[rno].getNewRegionNo();
					if(seedNo >0)
						regMatrix.set(i, j,seedNo);
				}
				int vno= (int)vMap.get(i, j);
				int seedNo = objProp[vno].getNewRegionNo();
				if(seedNo >0)
					vorMatrix.set(i, j,seedNo);
			}
		vorMatrix.removeHoles(Utility.FOUR);
		return vorMatrix;*/


	}



	private static Matrix regionGrowing(ImageObject[] objProp, Matrix objects, double thr) {

		int i,j,k;
		Matrix simMatrix = new Matrix(objProp.length,objProp.length,0.0);
		Matrix outMatrix = new Matrix(objects.getRowDimension(),objects.getColumnDimension(),0.0);
		int maxMatrixEntry = objects.maxMatrixEntry();
		for(i=1 ; i<maxMatrixEntry+1;++i){	
			ImageObject obj = objProp[i];
			for(k=0 ; k<obj.noNeighbours;++k){
				if(obj.getSimilarity()[k]<thr){
					simMatrix.set(i, obj.getNeigh()[k], 1);
				}
			}
		}


		simMatrix=	findConnectedComponents(simMatrix);

		for (i = 0; i < objects.getRowDimension(); i++)
			for (j = 0; j < objects.getColumnDimension(); j++){
				int rno= (int)objects.get(i, j);
				if(rno>0){
					int newRno = (int)simMatrix.get(rno, 0);
					outMatrix.set(i, j, newRno);
				}
			}
		return outMatrix;

	}

	public static Matrix findConnectedComponents(Matrix adjMatrix){
		Matrix groupNO = new Matrix(adjMatrix.getRowDimension(), 1, -1);
		int curGroupNO = 1;

		//start from the first node
		for(int i=0; i<adjMatrix.getRowDimension(); ++i){
			//if the node is part of a group, skip it
			if(groupNO.get(i, 0) != -1)
				continue;

			groupNO.set(i, 0, curGroupNO);

			int[] nodes = new int[adjMatrix.getRowDimension()];
			nodes[0] = i;
			int nodes_start = 0, nodes_end = 0, cur = 0;

			for(nodes_start=0; nodes_start<=nodes_end; nodes_start++){
				cur = nodes[nodes_start];

				for(int k=0; k<adjMatrix.getColumnDimension(); k++){
					if(adjMatrix.get(cur, k) > 0 && 
							groupNO.get(k, 0) == -1){
						groupNO.set(k, 0, curGroupNO);
						nodes[++nodes_end] = k;
					}
				}
			}

			curGroupNO++;
		}

		return groupNO;
	}

	public static double[] getObjFeatureVector(Matrix[] runLengthMatrices, ImageObject[] objProp, double centX, double centY, double numOfNodes, int radius, int featNo){
		double[][] matrices = new double[7][];
		double x,y;
		double dist;

		for(int connType=1; connType<=6; connType++)
			matrices[connType] = new double[runLengthMatrices[connType].getColumnDimension()];

		double countObj=0;
		for(int nodeID=1; nodeID<numOfNodes+1; nodeID++){

			ImageObject obj = objProp[nodeID];
			x = obj.centroidX;
			y = obj.centroidY;
			dist = Math.sqrt(Math.pow(centX-x, 2) + Math.pow(centY-y, 2));



			// if the node is within the circle, add its values to run
			// length matrix.


			if(dist < radius){
				countObj++;
				for(int connType=1; connType<=6; connType++)
					for(int col=0; col<runLengthMatrices[connType].getColumnDimension(); col++)
						matrices[connType][col] += ( runLengthMatrices[connType].get(nodeID, col));
			}

		}

		for(int connType=1; connType<=6; connType++)
			for(int col=0; col<runLengthMatrices[connType].getColumnDimension(); col++)
				matrices[connType][col] /= countObj;

		return getRunLengthFeatures(matrices, featNo, objProp);
	}

	private static double[] getRunLengthFeatures(double[][] matrices, int featNo, ImageObject[] objProp){

		//double [] features = new double[22];
		double [] features = new double[featNo];
		int i,j,index=0;

		/* n_r is the sum of run length numbers in the matrix*/
		double n_r = 0;
		double n_r_i[] = new double [7];
		for( i= 1 ; i< 7 ; ++i){
			for( j= 0 ; j< matrices[i].length ; ++j)	{
				n_r+=matrices[i][j];
				n_r_i[i]+=matrices[i][j];
			}
		}


		double edges =0; 
		double edges_i[] = new double [7];

		for(i=1; i<objProp.length;++i){
			ImageObject obj = objProp[i];

			for(j=0; j<obj.getNoNeighbours();++j){
				edges_i[obj.getType()[j]]++;
				edges++;
			} 
		}

		edges = edges/2;
		for(j=0; j<7;++j)
			edges_i[j]=edges_i[j]/2;

		/* short path emphasis*/
		double sum_all =0;
		for(i= 1 ; i< 7 ; ++i){
			double sum =0;
			for(j= 0 ; j< matrices[i].length ; ++j){			
				sum+= ( matrices[i][j] / Math.pow(j+1, 2) );	
				sum_all+= ( matrices[i][j] / Math.pow(j+1, 2) );
			}
			if(sum!=0 && n_r_i[i]!=0)
				features[index++]=sum/n_r_i[i]; //for a single edge type
			else
				features[index++]=0; //for a single edge type
		}
		if(featNo==16 || featNo==23 ){
			if(sum_all!=0 && n_r!=0)
				features[index++]=sum_all/n_r; 
			else
				features[index++]=0; 
		}
		//for all edge types

		/* long path emphasis*/
		sum_all =0;
		for(i= 1 ; i< 7 ; ++i){
			double sum =0;
			for(j= 0 ; j< matrices[i].length ; ++j){			
				sum+= ( matrices[i][j] * Math.pow(j+1, 2) );	
				sum_all+= ( matrices[i][j] * Math.pow(j+1, 2) );
			}
			if(sum!=0 && n_r_i[i]!=0)
				features[index++]=sum/n_r_i[i]; //for a single edge type
			else
				features[index++]=0; //for a single edge type

		}
		if(featNo==16 || featNo==23 ){
			if(sum_all!=0 && n_r!=0)
				features[index++]=sum_all/n_r; 
			else
				features[index++]=0; 
		}

		/*edge type nonuniformity*/
		sum_all =0;
		for(i= 1 ; i< 7 ; ++i){
			double sum =0;
			for(j= 0 ; j< matrices[i].length ; ++j){			
				sum+= ( matrices[i][j] ) ;	
			}
			//features[index++] = Math.pow(sum,2) /n_r ; //for a single edge type , /n_r is outside paranthesis
			sum_all+= Math.pow(sum,2);
		}
		if(sum_all!=0 && n_r!=0)
			features[index++]=sum_all/n_r; //for all edge types
		else
			features[index++]=0; 



		/*path length nonuniformity*/
		int max_length = 0;
		for (int k = 1; k < 7; k++)
			if(max_length < matrices[k].length)
				max_length = matrices[k].length;

		sum_all =0;
		for(j= 0 ; j< max_length ; ++j){
			double sum =0;
			for(i= 1 ; i< 7 ; ++i){		
				if(j < matrices[i].length)
					sum+= ( matrices[i][j] ) ;	
			}
			sum_all+= Math.pow(sum,2);
		}
		if(sum_all!=0 && n_r!=0)
			features[index++]=sum_all/n_r; //for all edge types
		else
			features[index++]=0;

		if(featNo==23){


			sum_all =0;
			for(i= 1 ; i< 7 ; ++i){
				double sum =0;
				for(j= 0 ; j< matrices[i].length ; ++j){			
					sum+= ( matrices[i][j]  );	
					sum_all+= ( matrices[i][j]  );
				}
				if(sum!=0 && n_r_i[i]!=0)
					features[index++]=sum/edges_i[i]; //for a single edge type
				else
					features[index++]=0; //for a single edge type
			}
			if(sum_all!=0 && n_r!=0)
				features[index++]=sum_all/edges; 
			else
				features[index++]=0; 

		}



		return features;

	}


	private static double[][] normalize(double[][] features) {

		int i, j;
		double N;
		double avg, std;

		for (i = 0; i < features[1].length; i++){
			avg = std = 0.0;
			N = 0.0;
			for (j = 1; j < features.length; j++)
			{
				avg += features[j][i];
				N++;
			}
			avg /= N;

			for (j = 1; j < features.length; j++)
				std += Utility.SQUARE(features[j][i] - avg);
			std = Math.sqrt(std / (N - 1));

			for (j = 1; j < features.length; j++)
				if ( (std > Utility.ZERO))
					features[j][i] = (features[j][i] - avg) / std;
		}
		return features;
	}


	private static void findSimlarityDistance(ImageObject[] objProp, Matrix objects) {

		int maxMatrixEntry = objects.maxMatrixEntry();
		for(int i=1 ; i<maxMatrixEntry+1;++i){
			ImageObject obj = objProp[i];
			for(int k=0 ; k<obj.noNeighbours;++k){
				double s = computeEuclideanDistance(obj,objProp[obj.getNeigh()[k]]);
				obj.getSimilarity()[k]=s;

			}
		}

	}

	private static double computeEuclideanDistance(ImageObject obj, ImageObject obj2) {
		int k;
		double d = 0; 

		for (k = 0; k < obj.ptex.length; k++)
			d += Utility.SQUARE(obj.ptex[k] - obj2.ptex[k]);
		return Math.sqrt(d);
	}




}