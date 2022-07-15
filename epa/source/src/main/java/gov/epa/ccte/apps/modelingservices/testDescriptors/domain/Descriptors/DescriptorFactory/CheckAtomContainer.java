package gov.epa.ccte.apps.modelingservices.testDescriptors.domain.Descriptors.DescriptorFactory;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;

public class CheckAtomContainer {

	//public static Set<String> supported = new HashSet<String>(Arrays.asList("C", "H", "O", "N", "F", "Cl", "Br", "I", "S", "P", "Si","As", "Hg", "Sn"));

	public static Set<String> supported = new HashSet<String>(Arrays.asList("H","B","C","N","O","F","Al","Si","P","S","Cl","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Mo","Ag","Cd","In","Sn","Sb","Te","I","Gd","Pt","Au","Hg","Tl","Pb","Bi"));
	
	public static void checkAtomContainer(IAtomContainer m, boolean also3D) throws IOException {

		if(also3D)HaveBadElement(m);

		if (m.getAtomCount() == 0)
			throw new IOException("Number of atoms equals zero");
		if (m.getAtomCount() == 1)
			throw new IOException("Only one nonhydrogen atom");
		if (!HaveCarbon(m)) 
			throw new IOException("Molecule does not contain carbon");


		AtomContainerSet  moleculeSet = (AtomContainerSet)ConnectivityChecker.partitionIntoMolecules(m);
		if (moleculeSet.getAtomContainerCount() > 1)
			throw new IOException("Multiple molecules");
	}

	static private void HaveBadElement(IAtomContainer mol) throws IOException {
		for (int i=0; i<mol.getAtomCount();i++) {
			String atom = mol.getAtom(i).getSymbol();
			if(!supported.contains(atom))
				throw new IOException("Molecule contains element: \"" + atom + "\" which is not supported for 3D WHIM descriptors");
		}

	}

	static private boolean HaveCarbon(IAtomContainer mol) {

		for (int i=0; i<mol.getAtomCount();i++) {
			String var = mol.getAtom(i).getSymbol();
			if (var.equals("C"))
				return true;
		}

		return false;

	}
}
