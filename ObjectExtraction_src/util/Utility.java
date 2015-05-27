package util;

public class Utility {

	public static final int FOUR = 4;
    public static final int EIGHT = 8;
    
    public static final int NODESHARE = 1;
    public static final int EDGESHARE = 2;
    
    public static final int UNVISITED = 0;
    public static final int VISITED = 1;
    public static final int NOT_IN = 2;
    
    public static final int LM_PIX = 3;
    public static final int NC_PIX = 1;
    public static final int ST_PIX = 2;
    
    public static final int BACKGROUND = 0;
    public static final int MYBACKGROUND = 0;
    public static final int UNSEGMENTED = 0;
    public static final int MARKED_UNSEGMENTED = -1;
    public static final int TEMPORARY_MARKED = -2;
    public static final int TEMP_VALUE = -400;
    public static final int NO_OBJECT = -1;
    
    public static final int LUMEN = 1;
    public static final int MAX_DISTANCE = 0;
    public static final int TOTAL_DISTANCE = 1;
    
    public static final double ZERO = 0.000000000000000001;
    public static final double SMALL = 0.00001;
    public static final double BIG = 9999999999.9;
    public static final double PI = 3.14159265358979;
    
    public static final int S_SQUARE = 0;
    public static final int S_OCTAGON = 1;
    public static final int CONSIDER_AREAS= 0;
    public static final int DO_NOT_CONSIDER_AREAS = 1;
	
	public static final int NUM_OF_CONN_TYPES = 6;
	
	public static final int EUCLIDEAN = 1;
	
	
	public static final int MAXDISTANCE = 4;

	public static final int ALLMERGE = -1;
	public static final int NOSEEDS = -2;
	public static final int NORMAL = 1;
	
	public static final int CANCER = 1;
	public static final int NONCANCER = 2;
	public static final int NOMATTER = 0;
	public static final int HE = 1;
	
	
    protected Utility(){
		//no body
	}
    
    public static int SQUARE(int a){
    	return a*a;
    }
    public static double SQUARE(double a){
    	return a*a;
    }
    
    public static int SQRDIST(int a, int b){
    	return ((a-b)*(a-b));
    }
    public static double SQRDIST(double a, double b){
    	return ((a-b)*(a-b));
    }
    
    public static void SWAP(int a, int b){
    	int temp=(a); (a) = (b); (b) = temp;
    }
    public static void SWAP(double a, double b){
    	double temp=(a); (a) = (b); (b) = temp;
    }

    public static void ROTATE(double a[][],int i,int j,int k,int l, double tau, double s){
    	double g=a[i][j]; 
    	double h=a[k][l]; 
    	a[i][j]=g-s*(h+g*tau); 
    	a[k][l]=h+s*(g-h*tau);
    }
    public static double hypot(double a, double b) {
        double r;
        if (Math.abs(a) > Math.abs(b)) {
           r = b/a;
           r = Math.abs(a)*Math.sqrt(1+r*r);
        } else if (b != 0) {
           r = a/b;
           r = Math.abs(b)*Math.sqrt(1+r*r);
        } else {
           r = 0.0;
        }
        return r;
     }
    
    public static int indexOfMinArrayEntry(double A[], int size){
    	int minIndex = 0, i;
    	for (i = 1; i < size; i++)
    		if (A[i] < A[minIndex])
    			minIndex = i;
    	return minIndex;
    }
    
    public static void jacobi(double a[][], int n, double d[], double v[][], int nrot[]){
        int j, iq, ip, i;
        double tresh, theta, tau, t, sm, s, h, g, c, b[], z[];
        
        b = new double[n];
        z = new double[n];
        for (ip = 0; ip < n; ip++){
            for (iq = 0; iq < n; iq++)
                v[ip][iq] = 0.0;
            v[ip][ip] = 1.0;
        }
        for (ip = 0; ip < n; ip++){
            b[ip] = d[ip] = a[ip][ip];
            z[ip]=0.0;
        }
        nrot[0] = 0;
        for (i = 1; i <= 50; i++){
            //printf("iteration %d\n",i);
            sm = 0.0;
            for (ip = 0; ip < n - 1; ip++){
                for (iq = ip + 1; iq < n; iq++)
                    sm += Math.abs(a[ip][iq]);
            }
            if (sm == 0.0){
                return;
            }
            if (i < 4)
                tresh = 0.2 * sm / (n * n);
            else
                tresh = 0.0;
            for (ip = 0; ip < n - 1; ip++){
                for (iq = ip + 1; iq < n; iq++){
                    g = 100.0 * Math.abs(a[ip][iq]);
                    if (i > 4 && (double)(Math.abs(d[ip]) + g) == (double)Math.abs(d[ip]) && 
                        (double)(Math.abs(d[iq]) + g) == (double)Math.abs(d[iq]))
                        a[ip][iq] = 0.0;
                    else if (Math.abs(a[ip][iq]) > tresh){
                        h = d[iq] - d[ip];
                        if ((double)(Math.abs(h) + g) == (double)Math.abs(h))
                            t = (a[ip][iq]) / h;
                        else{
                            theta = 0.5 * h / (a[ip][iq]);
                            t = 1.0 / (Math.abs(theta) + Math.sqrt(1.0 + theta * theta));
                            if (theta < 0.0)
                                t = -t;
                        }
                        c = 1.0 / Math.sqrt(1 + t * t);
                        s = t * c;
                        tau = s / (1.0 + c);
                        h = t * a[ip][iq];
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= h;
                        d[iq] += h;
                        a[ip][iq] = 0.0;
                        for (j = 0; j <= ip - 1; j++){
                            ROTATE(a,j,ip,j,iq,tau,s);
                        }
                        for (j = ip + 1; j <= iq - 1; j++){
                            ROTATE(a,ip,j,j,iq,tau,s);
                        }
                        for (j = iq + 1; j < n; j++){
                            ROTATE(a,ip,j,iq,j,tau,s);
                        }
                        for (j = 0; j < n; j++){
                            ROTATE(v,j,ip,j,iq,tau,s);
                        }
                        ++(nrot[0]);
                    }
                }
            }
            for (ip = 0; ip < n; ip++){
                b[ip] += z[ip];
                d[ip] = b[ip];
                z[ip] = 0.0;
            }       
        }
        System.out.println("too many iterations in routine jacobi");
    }
    public static void eigsrt(double d[], double v[][], int n){
        int k, j, i;
        double p;
        
        for (i = 0; i < n; i++){
            p = d[k = i];
            for (j = i + 1; j < n; j++)
                if (d[j] >= p)
                    p = d[k = j];
            if (k != i){
                d[k] = d[i];
                d[i] = p;
                for (j = 0; j < n; j++){
                    p = v[j][i];
                    v[j][i] = v[j][k];
                    v[j][k] = p;
                }
            }
        }
    }
    
    public static void gaussj(double a[][], int n, double b[][], int m){
    	int indxc[], indxr[], ipiv[];
    	int i, icol=0, irow=0, j, k, l, ll;
    	double big, dum, pivinv;

    	indxc = new int[n];
    	indxr = new int[n];
    	ipiv = new int[n];
    	for (j = 0; j < n; j++)
    		ipiv[j] = 0;
    	for (i = 0; i < n; i++){
    		big = 0.0;
    		for (j = 0; j < n; j++)
    			if (ipiv[j] != 1)
    				for (k = 0; k < n; k++){
    					if (ipiv[k] == 0){
    						if (Math.abs(a[j][k]) >= big){
    							big = Math.abs(a[j][k]);
    							irow = j;
    							icol = k;
    						}
    					}
    					else if (ipiv[k] > 1){
    						System.out.println("\nError: gaussj: Singular matrix - 1\n");
    						System.exit(1);
    					}
    				}
    		++(ipiv[icol]);
    		if (irow != icol){
    			for (l = 0; l < n; l++)
    				SWAP(a[irow][l],a[icol][l]);
    			for (l = 0; l < m; l++)
    				SWAP(b[irow][l],b[icol][l]);
    		}
    		indxr[i] = irow;
    		indxc[i] = icol;
    		if (a[icol][icol] == 0.0){
    			System.out.println("\nError: gaussj: Singular matrix - 2\n");
    			System.exit(1);
    		}
    		pivinv = 1.0 / a[icol][icol];
    		a[icol][icol] = 1.0;
    		for (l = 0; l < n; l++)
    			a[icol][l] *= pivinv;
    		for (l = 0; l < m; l++)
    			b[icol][l] *= pivinv;
    		for (ll = 0; ll < n; ll++)
    			if (ll != icol){
    				dum = a[ll][icol];
    				a[ll][icol] = 0.0;
    				for (l = 0; l < n; l++)
    					a[ll][l] -= a[icol][l] * dum;
    				for (l = 0; l < m; l++)
    					b[ll][l] -= b[icol][l] * dum;
    			}
    	}
    	for (l = n - 1; l >= 0; l--){
    		if (indxr[l] != indxc[l])
    			for (k = 0; k < n; k++)
    				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    	}
    	
    }
}








