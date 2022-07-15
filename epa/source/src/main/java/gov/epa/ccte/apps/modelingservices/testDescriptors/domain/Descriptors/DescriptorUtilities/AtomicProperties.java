package gov.epa.ccte.apps.modelingservices.testDescriptors.domain.Descriptors.DescriptorUtilities;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Hashtable;
import java.util.LinkedList;

import qspr.metaserver.cs.EPAServer;


public class AtomicProperties {

	private static AtomicProperties ap=null;
	
	private Hashtable htMass=new Hashtable();
	private Hashtable htVdWVolume=new Hashtable();
	private Hashtable htElectronegativity=new Hashtable();
	private Hashtable htPolarizability=new Hashtable();
	
	
	private AtomicProperties() throws IOException {
						
		BufferedReader br=new BufferedReader(new FileReader(EPAServer.getFile("whim_weights.txt")));
		
		String Header=br.readLine(); // header
		
		String Line="";
		while (true) {
			Line=br.readLine();
			if (!(Line instanceof String)) {
				break;
			}
			
			LinkedList l=Utilities.Parse(Line,"\t");
			
			String symbol=(String)l.get(0);
			htMass.put(symbol,l.get(1));
			htVdWVolume.put(symbol,l.get(2));
			htElectronegativity.put(symbol,l.get(3));
			htPolarizability.put(symbol,l.get(4));
			
		}
						
		br.close();
		
	}

	public double GetVdWVolume(String symbol) {
		double VdWVolume=-99;
		
		String strVdWVolume=(String)htVdWVolume.get(symbol);
		
		try {
			VdWVolume=Double.parseDouble(strVdWVolume);
		} catch (Exception e) {
			System.out.println(e);
		}
		
		
		return VdWVolume;
		
	}
	
	public double GetNormalizedVdWVolume(String symbol) {
		double VdWVolume=-99;
		
		VdWVolume=this.GetVdWVolume(symbol)/this.GetVdWVolume("C");
				
		return VdWVolume;
		
	}
	
	public double GetElectronegativity(String symbol) {
		double Electronegativity=-99;
		
		String strElectronegativity=(String)htElectronegativity.get(symbol);
		
		try {
		Electronegativity=Double.parseDouble(strElectronegativity);
		} catch (Exception e) {
			System.out.println(e);
		}
		
		
		return Electronegativity;
		
	}
	
	public double GetNormalizedElectronegativity(String symbol) {
		double Electronegativity=-99;
		
		Electronegativity=this.GetElectronegativity(symbol)/this.GetElectronegativity("C");
				
		return Electronegativity;
		
	}
	public double GetPolarizability(String symbol) {
		double Polarizability=-99;
		
		String strPolarizability=(String)htPolarizability.get(symbol);
		
		try {
		Polarizability=Double.parseDouble(strPolarizability);
		} catch (Exception e) {
			System.out.println(e);
		}
		
		
		return Polarizability;
		
	}
	
	public double GetNormalizedPolarizability(String symbol) {
		double Polarizability=-99;
		
		Polarizability=this.GetPolarizability(symbol)/this.GetPolarizability("C");
				
		return Polarizability;
		
	}
	public double GetMass(String symbol) {
		double mass=-99;
		
		String strMass=(String)htMass.get(symbol);
		
		try {
		mass=Double.parseDouble(strMass);
		
		} catch (Exception e) {
			System.out.println("AtomicProperties--GetMass:"+symbol);
		}
		
		
		return mass;
		
	}
	
	public double GetNormalizedMass(String symbol) {
		double mass=-99;
		
		mass=this.GetMass(symbol)/this.GetMass("C");
				
		return mass;
		
	}
	
	
	
	public static AtomicProperties getInstance() throws IOException
	{
		if (ap == null) {
			ap = new AtomicProperties();
		}
		return ap;
	}
}
