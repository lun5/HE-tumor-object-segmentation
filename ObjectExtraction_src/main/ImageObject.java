package main;


import java.io.BufferedReader;
import java.util.StringTokenizer;

import matrix.*;
import util.*;

public class ImageObject {

	double area;
	double centroidX, centroidY;
	int regionNo;
	int objectID=0;
	

	int newRegionNo;
	int mergeRegionNo;
	int clusterNo;
	int superNode;
	double[] ptex;
	int noNeighbours;
	int[] neigh;
	int[] type;
	double[] similarity;
	int index;
	boolean isBorder;
	boolean isSelected;
	int cancer;
	int noTriangles;
	public ImageObject(int regionNo){
		this.area = 0.0;
		this.centroidX = 0.0;
		this.centroidY = 0.0;
		this.regionNo = regionNo;
		this.newRegionNo = -1;
		this.noNeighbours = 0;
		this.index=0;
		this.objectID=regionNo;
		this.superNode=0;
		noTriangles=0;
		isBorder= false;
		isSelected=false;
	}
	
	public ImageObject(double area, double centroidX, double centroidY, int regionNo) {
		this.area = area;
		this.centroidX = centroidX;
		this.centroidY = centroidY;
		this.regionNo = regionNo;
		this.objectID=regionNo;
		this.noNeighbours = 0;
		this.superNode=0;
		isBorder= false;
	}



	public static ImageObject[] computeObjectPropertiesL(Matrix cmap, Matrix objects, int objectNo) {
		// TODO Auto-generated method stub
//		 regions start with label 1
		ImageObject objProp[] = new ImageObject[objectNo+1];
		int i, j, rno, cno;
		for (i = 1; i <= objectNo; i++){
			objProp[i]= new ImageObject(i);
			
		}

		for (i = 0; i < objects.getRowDimension(); i++)
			for (j = 0; j < objects.getColumnDimension(); j++){
				/*if (objects.get(i, j) == Utility.BACKGROUND){
					System.out.println("\n\nError: BACKGROUND in computeObjectProperties\n\n");
					System.exit(1);
				}*/
				rno = (int)objects.get(i, j);
				if(rno==Utility.BACKGROUND) continue;
				
				cno = (int)cmap.get(i, j);

				if(cno==Utility.MYBACKGROUND) continue;
				
				objProp[rno].area++;
				objProp[rno].centroidX += i;
				objProp[rno].centroidY += j;
				objProp[rno].clusterNo = cno;
			}
		for (i = 1; i <= objectNo; i++){
			if (objProp[i].area > 0){
				objProp[i].centroidX /= objProp[i].area;
				objProp[i].centroidY /= objProp[i].area;
			}
			
		}
		return objProp;
		
	}

	public double getArea() {
		return area;
	}

	public void setArea(double area) {
		this.area = area;
	}

	public double getCentroidX() {
		return centroidX;
	}

	public void setCentroidX(double centroidX) {
		this.centroidX = centroidX;
	}

	public double getCentroidY() {
		return centroidY;
	}

	public void setCentroidY(double centroidY) {
		this.centroidY = centroidY;
	}

	public int getClusterNo() {
		return clusterNo;
	}

	public void setClusterNo(int clusterNo) {
		this.clusterNo = clusterNo;
	}

	public int getRegionNo() {
		return regionNo;
	}

	public void setRegionNo(int regionNo) {
		this.regionNo = regionNo;
	}

	public int getNewRegionNo() {
		return newRegionNo;
	}

	public void setNewRegionNo(int newRegionNo) {
		this.newRegionNo = newRegionNo;
	}
	public void createNeighs(int neighNo) {
		neigh= new int[neighNo];
		type= new int[neighNo];
		similarity = new double[neighNo];
	}

	public static void readList(ImageObject[] objProp, BufferedReader input) throws java.io.IOException{
		StringTokenizer ST;
		String str;
		int n=0;
		if ((str=input.readLine())!=null){
			ST = new StringTokenizer(str);
			n=(int)Double.parseDouble(ST.nextToken());
		}
		for(int i=0; i<n;++i){
			 
			str=input.readLine();
			ST = new StringTokenizer(str," ");
			int index= Integer.parseInt(ST.nextToken());
			int neighNo= Integer.parseInt(ST.nextToken());
			objProp[index].createNeighs(neighNo);
			objProp[index].noNeighbours=neighNo;
			for(int j=0; j<neighNo; ++j){
				objProp[index].neigh[j]=Integer.parseInt(ST.nextToken());
				objProp[index].type[j]=Integer.parseInt(ST.nextToken());
			}
		}
		
	}

	public int[] getNeigh() {
		return neigh;
	}

	public void setNeigh(int[] neigh) {
		this.neigh = neigh;
	}

	public int getNoNeighbours() {
		return noNeighbours;
	}

	public void setNoNeighbours(int noNeighbours) {
		this.noNeighbours = noNeighbours;
	}

	public double[] getPtex() {
		return ptex;
	}

	public void setPtex(double[] ptex) {
		this.ptex = ptex;
	}

	public int[] getType() {
		return type;
	}

	public void setType(int[] type) {
		this.type = type;
	}

	public double[] getSimilarity() {
		return similarity;
	}

	public void setSimilarity(double[] similarity) {
		this.similarity = similarity;
	}

	public void updateNeighs(ImageObject object) {

		if(noNeighbours==0){
			neigh= new int[++noNeighbours];
			type= new int[noNeighbours];
			this.neigh[this.index]=object.getRegionNo();
			this.type[this.index]=getEdgeType(object);
			similarity = new double[noNeighbours];
		}
		else{
			int control = 0;
			for(int i=0 ; i<neigh.length; ++i){
				if(neigh[i]==object.getRegionNo()){
					control=1;
					break;
				}		
			}
			if(control!=1){
				int[] neigh2= new int[noNeighbours];
				int[] type2= new int[noNeighbours];
				for(int i=0 ; i<neigh.length; ++i){
					neigh2[i]=neigh[i];
					type2[i]=type[i];
				}
				neigh= new int[++noNeighbours];
				type= new int[noNeighbours];
				similarity = new double[noNeighbours];
				for(int i=0 ; i<neigh2.length; ++i){
					neigh[i]=neigh2[i];
					type[i]=type2[i];
				}
				this.neigh[++this.index]=object.getRegionNo();
				this.type[this.index]=getEdgeType(object);
			}
		}
	}
	public int getObjectID() {
		return objectID;
	}

	public void setObjectID(int objectID) {
		this.objectID = objectID;
	}
	public int getEdgeType(ImageObject obj){
		if(obj.getClusterNo()==0 && this.getClusterNo()==0){
			return 1;
		}
		else if(obj.getClusterNo()==1 && this.getClusterNo()==1){
			return 2;
		}
		else if(obj.getClusterNo()==2 && this.getClusterNo()==2){
			return 3;
		}
		else if((obj.getClusterNo()==1 && this.getClusterNo()==0) || 
				(obj.getClusterNo()==0 && this.getClusterNo()==1) ){
			return 4;
		}
		else if((obj.getClusterNo()==0 && this.getClusterNo()==2) || 
				(obj.getClusterNo()==2 && this.getClusterNo()==0) ){
			return 5;
		}
		else if((obj.getClusterNo()==1 && this.getClusterNo()==2) || 
				(obj.getClusterNo()==2 && this.getClusterNo()==1) ){
			return 6;
		}
		return 0;
		
	}

	public int getMergeRegionNo() {
		return mergeRegionNo;
	}

	public void setMergeRegionNo(int mergeRegionNo) {
		this.mergeRegionNo = mergeRegionNo;
	}


	public int getSuperNode() {
		return superNode;
	}

	public void setSuperNode(int superNode) {
		this.superNode = superNode;
	}

	public boolean isSelected() {
		return isSelected;
	}

	public void setSelected(boolean isSelected) {
		this.isSelected = isSelected;
	}
}
