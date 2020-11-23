//Servidor
import java.io.*;
import java.net.*;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.StringTokenizer;
import java.util.UUID;
import java.util.logging.Level;
import java.util.logging.Logger;

class CommandSlaver {

	static class ConnectionHandler implements Runnable {
		private Socket clsSocket;
		private DataInputStream clsInputStream;
		private DataOutputStream clsOutputStream;		

		private Spitter clsSpitter;
		private Thread clsSpitterThread;
		
		private Process clsProcess;	
		private OutputStreamWriter ow;
		private InputStreamReader isr;
		
		private String clsBench;
		
		private UUID clsSession;
		private String clsSessionPath;
		
		static class Spitter implements Runnable {
			private DataOutputStream clsOutputStream;	
			private InputStream clsInputStream;
						
			Spitter(InputStream i, DataOutputStream o) {
				clsOutputStream = o;
				clsInputStream = i;
			} // end constructor
			
			public void run() {
				try {
					BufferedReader br = new BufferedReader(new InputStreamReader(clsInputStream));
					
					while (true) {
						String line = br.readLine();
						System.out.println(line);

						clsOutputStream.writeBytes(line);
						clsOutputStream.writeByte('\n');
					} // end while
					//clsOutputStream.writeChar('\0');					
				} catch( Exception ex ){
					Logger.getLogger(CommandSlaver.class.getName()).log(Level.SEVERE, null, ex);
				} // end catch
			} // end method
		} // end class
		
		public ConnectionHandler(Socket socket) throws IOException {
			clsSocket = socket;
			clsInputStream = new DataInputStream(clsSocket.getInputStream());
			clsOutputStream = new DataOutputStream(clsSocket.getOutputStream());
			
			
			
		} // end constructor
		
		public void run() {
			try {
				
				clsSession = UUID.randomUUID();
				clsSessionPath = clsSession + "/";
				
				createDirectory(clsSessionPath);
				copyFile( "stuff/run.sh", clsSessionPath);
				copyFile( "stuff/pt_load_scripts.tcl", clsSessionPath);
				
				File script = new File(clsSessionPath + "run.sh");
				script.setExecutable(true);
								
				while (true) {
					StringBuffer command = new StringBuffer();
					
					System.out.print("Waiting for new command line...");
					
					int ch;
					while ( (ch = clsInputStream.read()) != '\n') {
						System.out.print((char)ch);
						command.append((char) ch);
					} // end while

					System.out.print("Command Line: " + command + '\n');

					if ( !parseCommand(command.toString()) )
						break;
				} // end while
				
				System.out.print("Connection CLOSED");
				
				// Close connection.
				clsInputStream.close();
				clsOutputStream.close();

			} catch (IOException ex) {
				Logger.getLogger(CommandSlaver.class.getName()).log(Level.SEVERE, null, ex);
			} catch (Exception ex ) {
				Logger.getLogger(CommandSlaver.class.getName()).log(Level.SEVERE, null, ex);
			} // end catch
		} // end method
		
		// ---------------------------------------------------------------------
		
		@Override
		protected void finalize() throws Throwable {
			super.finalize();
			
			deleteDirectory(new File(clsSessionPath));
		} // end method
		
		// ---------------------------------------------------------------------
		
		private void createDirectory(String path) throws IOException {
			try {
				Files.createDirectory(Paths.get(path));
			} catch ( FileAlreadyExistsException e ) {
			} // end catch
		} // end method		

		// ---------------------------------------------------------------------
		
		// Source: http://stackoverflow.com/questions/7768071/java-delete-a-folder-content
		
		private void deleteDirectory(File folder) throws IOException {
			File[] files = folder.listFiles();
			if (files != null) {
				for (File f : files) {
					if (f.isDirectory()) {
						deleteDirectory(f);
					} else {
						f.delete();
					} // end else
				} // end for
			} // end if
			folder.delete();
		} // end method
		
		// ---------------------------------------------------------------------

		private void copyFile(String source, String target) throws IOException {
			// [NOTE] Using getResourceAsStream() and not a simple file copy since
			// in the realease version the source file will be insite the jar.

			if ( Files.isDirectory(Paths.get(target)) )
				target += "/" + Paths.get(source).getFileName();
			
			InputStream inputStream = this.getClass().getResourceAsStream(source);
			OutputStream out = new FileOutputStream(new File(target));

			int read = 0;
			byte[] bytes = new byte[1024];

			while ((read = inputStream.read(bytes)) != -1) {
				out.write(bytes, 0, read);
			} // end while

			inputStream.close();
			out.flush();
			out.close();
		} // end method		
		
		// ---------------------------------------------------------------------
		
		private void performStart(StringTokenizer tokens) throws IOException {
			clsBench = tokens.nextToken();
			clsProcess = Runtime.getRuntime().exec("./run.sh " + clsBench, null, new File(clsSessionPath));

			ow = new OutputStreamWriter(clsProcess.getOutputStream());
				
			clsSpitter = new Spitter(clsProcess.getInputStream(), clsOutputStream);
			clsSpitterThread = new Thread(clsSpitter);
			clsSpitterThread.start();
		} // end method
		
		// ---------------------------------------------------------------------
		
		private void performClose(StringTokenizer tokens) throws IOException {
			clsProcess.destroy();
		} // end method		
		
		// ---------------------------------------------------------------------
		
		private void performExec(StringTokenizer tokens) {
			try {

				StringBuffer command = new StringBuffer();

				int ch;
				while ((ch = clsInputStream.read()) != '\n') {
					System.out.print((char) ch);
					command.append((char) ch);
				} // end while

				ow.write(command.toString() + "\n");
				ow.flush();
			} catch (Exception e) {
				sendMessage("Exception: " + e.getMessage());
			} // end catch
		} // end method

		// ---------------------------------------------------------------------

		private void performDownload(StringTokenizer tokens) throws Exception {
			byte[] b = new byte[1024];
			int len = 0;
			int byteCount =0;
			String FileName = tokens.nextToken();

			FileInputStream file = new FileInputStream(FileName);

			while ( (len = file.read(b, 0, 1024)) != -1) {
				byteCount += len;
				clsOutputStream.write( b, 0, len);
			} // end while
			file.close();
			System.out.println("\tFile Sent : " + byteCount + " bytes" );		
		} // end method

		// ---------------------------------------------------------------------

		private void performUpload(StringTokenizer tokens) throws Exception {

			int len = 0;
			int byteCount =0;

			int length = Integer.parseInt(tokens.nextToken());

			FileOutputStream file = new FileOutputStream(clsSessionPath + clsBench + ".int.sizes");

			for ( int i = 0; i < length; i++ ) {
				file.write(clsInputStream.read());
			} // end for
			//System.out.print(new String(b));

			file.close();
			System.out.println("File Received : " + length + " bytes" );		
		} // end method

		// ---------------------------------------------------------------------

		private boolean parseCommand(String command) throws Exception {
			StringTokenizer tokens = new StringTokenizer(command);
			String operation = tokens.nextToken();

			if (operation.compareTo("start") == 0) {
				performStart(tokens);
				return true;
			} else if (operation.compareTo("close") == 0) {
				return false;		
			} else if (operation.compareTo("exec") == 0) {
				performExec(tokens);
				return true;
			} else if (operation.compareTo("download") == 0) {
				performDownload(tokens);
				return true;
			} else if (operation.compareTo("upload") == 0) {
				performUpload(tokens);
				return true;
			} else {
				sendMessage( "Invalid command.");
				return true;
			} // end else
		} // end method

		// ---------------------------------------------------------------------

		private void sendMessage( String msg ) {
			try {
				clsOutputStream.writeBytes(msg);
				clsOutputStream.writeByte( '\0');
			} catch ( Exception e ) {
				System.err.println( "Exception: "  + e.getMessage() );
			} // end catch
		} // end method
		
	} // end class
	
	public static void main(String argv[]) throws Exception {
		if ( argv.length != 1 ) {
			System.err.println( "Usage: java CommandServer <port number>");
			System.exit(1);
		} // end if
		
		int portNumber = Integer.parseInt( argv[0] );
		
		// creat a socket
		ServerSocket welcomeSocket = new ServerSocket( portNumber );

		// wait for clients
		while (true) {
			try {

				// wait for a connection
				Socket socket = welcomeSocket.accept();
				//System.out.println( "Connection accepted!" );	    

				ConnectionHandler handler = new ConnectionHandler(socket);
				
				Thread thread = new Thread(handler);
				thread.start();
			} catch ( Exception e ) {
				System.err.println( "Exception: " + e.getMessage() );
			} // end catch
		} // end while
	} // end main


	
} // end class
