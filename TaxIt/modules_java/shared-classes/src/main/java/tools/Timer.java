/* Copyright (c) 2016,
 * Mathias Kuhring, KuhringM@rki.de, Robert Koch Institute, Germany,
 * All rights reserved. For details, please note the license.txt.
 */

package tools;

public class Timer {
	// simple timer, use with timerStart() and timerStop()
	public long timer = 0;

	// starts a simple timer
	public void start(){
		timer = System.currentTimeMillis();
	}

	// stops the simple timer
	public void stop(){
		if (timer == 0)  System.out.println("warning: timer not startet");
		else System.out.print((System.currentTimeMillis() -  timer)/1000 + " sec");	
	}
}
