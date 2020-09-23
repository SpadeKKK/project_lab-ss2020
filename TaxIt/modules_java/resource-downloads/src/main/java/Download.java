import java.io.BufferedReader;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.SocketException;
import java.net.URL;

public class Download {

	// test strings for old HTTP response code handling
//	private static final List<String> IOE_TESTS = Arrays.asList(
//			"Server returned HTTP response code: 404 for URL:",		// java.io.FileNotFoundException ?
//			"Server returned HTTP response code: 502 for URL:",
//			"Server returned HTTP response code: 503 for URL:");

	// first not acceptable HTTP response code
	private static final short MAX_CODE = 400;

	// waiting time for repeating after a bad connection (e.g. HTTP response code >= MAX_CODE)
	private static final int WAIT_MIN = 1;
	private static final int WAIT_MILL = WAIT_MIN * 60 * 1000;

	// size of download pieces for the BufferedReader
	private static final int BUFFER_SIZE=1024;

	// changed response code handling based on:
	// http://stackoverflow.com/a/17352452
	// http://stackoverflow.com/a/4136963

	// which response codes should lead to a repeat, which to an exit?
	// https://de.wikipedia.org/wiki/HTTP-Statuscode
	
	// Similar problems with NCBI Entrez in Python: 
	// https://www.biostars.org/p/121287/

	// download data/file from a given url (String), return as String
	// queries will be repeated in case of HTTP response codes >= 400
	// (which may occur occasionally and usually temporarily in case of frequent queries to e.g. NCBI databases)
	public static String downloadFile(String url){
		StringBuffer sb = null;
		boolean done = false;

		// repeat until download was succesful (or an unhandled IOException occurs)
		while(!done){
			sb =  new StringBuffer();

			try {
				// establish connection
				HttpURLConnection httpConn = (HttpURLConnection) new URL(url).openConnection();

				// check for bad HTTP response codes
				if (httpConn.getResponseCode() >= MAX_CODE) {
					System.err.println("error: server problem (HTTP response code: " 
				    		+ httpConn.getResponseCode() 
				    		+  "), retry after " + WAIT_MIN + " min...");
					System.err.println(url);
					waitBeforeNextTry(WAIT_MILL);
				} 
				else {
					// pass connection to reader
					BufferedReader br = new BufferedReader(new InputStreamReader(httpConn.getInputStream()));

					// piecewise download based on http://stackoverflow.com/a/13210179
					char[] buffer = new char[BUFFER_SIZE];
					int charsRead = 0;
					try {
						while ((charsRead = br.read(buffer, 0, BUFFER_SIZE)) != -1) {
							sb.append(buffer, 0, charsRead);
						}
					}
					finally{
						br.close();
						done = true;
					}
				}
			}
			catch (EOFException e) {
				// java.io.IOException: Premature EOF -> seems not to be captured by EOFException but IOException
				System.err.println("error: download problem (EOFException), retry after " + WAIT_MIN + " min...");
				System.err.println(url);
				System.err.println(e);
				waitBeforeNextTry(WAIT_MILL);
			}
			catch (SocketException e) {
				System.err.println("error: download problem (SocketException), retry after " + WAIT_MIN + " min...");
				System.err.println(url);
				System.err.println(e);
				waitBeforeNextTry(WAIT_MILL);
			}
			catch (IOException e) {
				// hoped to catch "Premature EOF" with EOFException or this: http://stackoverflow.com/a/13210179
				// didn't work out, so now it is catched via IOException and string matching
				if (e.getMessage().startsWith("Premature EOF")){
					System.err.println("error: download problem (IOException: Premature EOF), retry after " + WAIT_MIN + " min...");
					System.err.println(url);
					System.err.println(e);
					waitBeforeNextTry(WAIT_MILL);
				}
				else {
					e.printStackTrace();
					System.err.println("error: download problem (unknown IOException), retry after " + WAIT_MIN + " min...");
					System.err.println(url);
					System.err.println(e);
					waitBeforeNextTry(WAIT_MILL);
				}

				// old HTTP response code catching based on string matching on the exception message
//				boolean repeat = false;
//				for (String test : IOE_TESTS){
//					repeat = repeat || e.getMessage().startsWith(test);
//				}
//				if (repeat){
//					System.err.println("error: server problem, retry...");
//				}
//				else {
//					e.printStackTrace();
//					System.err.println("error: downloads failed");
//					System.err.println(e.getMessage());
//					System.exit(1);
//				}
			}
		}
		return sb.toString();
	}

	// give the server a break before repeating a download
	private static void waitBeforeNextTry(int millSecs){
		try {
			Thread.sleep(millSecs); }
		catch (InterruptedException ie) {
			ie.printStackTrace();
		}
	}
}
