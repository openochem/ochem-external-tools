package gov.epa.ccte.apps.modelingservices.testDescriptors.domain.Descriptors.DescriptorUtilities;

import org.openscience.cdk.interfaces.IAtomContainer;


public class Inchi {

	public String inchi, inchiKey, inchiKey1,warning; 
	
	public static Inchi generateInChiKeyCDK(IAtomContainer ac) {		
		return CDKUtilities.generateInChiKey(ac);		
	}
	
	public static Inchi generateInChiKeyIndigo(IAtomContainer ac) {		
		return IndigoUtilities.generateInChiKey(ac);		
	}

	public static Inchi generateInChiKeyIndigo(String smiles) {
		return IndigoUtilities.toInchiIndigo(smiles);	
	}
	
	
}
