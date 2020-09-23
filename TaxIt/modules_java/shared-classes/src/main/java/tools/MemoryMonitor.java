/* Copyright (c) 2016,
 * Mathias Kuhring, KuhringM@rki.de, Robert Koch Institute, Germany,
 * All rights reserved. For details, please note the license.txt.
 */

package tools;

// simple memory monitor class
// source: http://viralpatel.net/blogs/getting-jvm-heap-size-used-memory-total-memory-using-java-runtime/
public class MemoryMonitor {
    
	private Runtime runtime;
	
	public MemoryMonitor(){
		//Getting the runtime reference from system
		runtime = Runtime.getRuntime();
	}

    // print memory usage in GB
    public void printGbUsage() {
         int scale = 1024*1024*1024;
         print(scale);
    }

    // print memory usage in bytes/scale
    private void print(int scale){      
        System.out.println("##### Heap utilization statistics [GB] #####");
         
        //Print used memory
        System.out.println("Used Memory:\t"
            + (runtime.totalMemory() - runtime.freeMemory()) / scale);
 
        //Print free memory
        System.out.println("Free Memory:\t"
            + runtime.freeMemory() / scale);
         
        //Print total available memory
        System.out.println("Total Memory:\t" + runtime.totalMemory() / scale);
 
        //Print maximum available memory
        System.out.println("Max Memory:\t" + runtime.maxMemory() / scale);
    }
}