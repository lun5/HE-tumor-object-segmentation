// HelloWorld.java

public class HelloWorld{

	public HelloWorld(){
		setVersion(0);
	}

	public static void main(String[] args){

		System.out.println("Hello World!");
	}

	public void setVersion(int aVersion){
		if( aVersion < 0 ){
			System.err.println("Improper version specified.");
		}
		else{
			version = aVersion;
		}
	}

	public int getVersion(){
		return version;
	}


	private int version;

}