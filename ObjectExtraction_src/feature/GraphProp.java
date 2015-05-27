package feature;

import main.ImageObject;
import matrix.*;
import util.*;
import java.util.Vector;

public class GraphProp {

	public static IntMatrix[] findConnectedComponents(IntMatrix adjMatrix){
		IntMatrix[] connectedComponents  = new IntMatrix[6];

		for(int connType=1; connType<=Utility.NUM_OF_CONN_TYPES; connType++){
			IntMatrix groupNO = new IntMatrix(adjMatrix.row, 1);
			groupNO.initializeMatrix(-1);
			connectedComponents[connType-1] = groupNO;
			int curGroupNO = 1;

			//start from the first node
			for(int i=0; i<adjMatrix.row; ++i){
				//if the node is part of a group, skip it
				if(groupNO.data[i][0] != -1)
					continue;

				groupNO.data[i][0] = curGroupNO;

				int[] nodes = new int[adjMatrix.row];
				nodes[0] = i;
				int nodes_start = 0, nodes_end = 0, cur = 0;

				for(nodes_start=0; nodes_start<=nodes_end; nodes_start++){
					cur = nodes[nodes_start];

					for(int k=0; k<adjMatrix.column; k++){
						if(adjMatrix.data[cur][k] == connType && 
								groupNO.data[k][0] == -1){
							groupNO.data[k][0] = curGroupNO;
							nodes[++nodes_end] = k;
						}
					}
				}

				curGroupNO++;
			}
		}

		return connectedComponents;
	}

	public static IntMatrix[] findConnCompLengths(IntMatrix adjMatrix){
		IntMatrix[] connComp = findConnectedComponents(adjMatrix);
		IntMatrix[] compLengths = new IntMatrix[6];

		for(int connType=1; connType<=Utility.NUM_OF_CONN_TYPES; connType++){
			IntMatrix curConnComp = connComp[connType-1];
			int numOfComp = curConnComp.maxElement();
			IntMatrix curCompLengths = new IntMatrix(numOfComp,1);

			//find the length of each component
			for(int i=0; i<curConnComp.row; i++)
				curCompLengths.data[(curConnComp.data[i][0]-1)][0] += 1;

			compLengths[connType-1] = curCompLengths;
		}

		return compLengths;
	}

	
	public static Matrix[] findMaxRunLengthMatrixForEachNode(ImageObject[] objProp, int winSize){
		Matrix[] runLengthMatrices = new Matrix[7];
		for(int i=1; i<=Utility.NUM_OF_CONN_TYPES; i++)
			runLengthMatrices[i] = new Matrix(objProp.length, 15);

		int[] distance, nodes;
		int nodes_start = 1, nodes_end = 1, cur = 1;

		//for each connection type...
		for(int connType=1; connType<=Utility.NUM_OF_CONN_TYPES; connType++){
			Matrix isUsed = new Matrix(objProp.length, 1, -1);

			//start from the first node
			for(int curRoot=1; curRoot<objProp.length; ++curRoot){
				nodes = new int[objProp.length];
				distance = new int[objProp.length];
				Vector<Integer> nodesInWindow = new  Vector<Integer>();
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
						if((int)objProp[cur].getType()[k] == connType && 
								(int)isUsed.get(objProp[cur].getNeigh()[k], 0) == -1 &&
								dist<winSize){
							//TODO inside window runlength control
							isUsed.set(objProp[cur].getNeigh()[k], 0, 0);
							
							nodesInWindow.add(objProp[cur].getNeigh()[k]);
							hasUnvisitedNeighbour = true; 

							nodes[++nodes_end] = objProp[cur].getNeigh()[k];
							distance[nodes_end] = distance[nodes_start] + 1;

						}
					}

					if(hasUnvisitedNeighbour == false){
						if(runLengthMatrices[connType].getColumnDimension() > distance[nodes_start]){
							runLengthMatrices[connType].set(curRoot, distance[nodes_start], 
									runLengthMatrices[connType].get(curRoot, distance[nodes_start])+1);

						}
						else{
							// enlarge the matrix by column if necessary
							runLengthMatrices[connType].enlargeByColumn(100);
							runLengthMatrices[connType].set(curRoot, distance[nodes_start], 
									runLengthMatrices[connType].get(curRoot, distance[nodes_start])+1);
						}
					}


				}

			}
		}

		return runLengthMatrices;
	}
	
	

}