package gov.epa.ccte.apps.modelingservices.testDescriptors.service;

import gov.epa.ccte.apps.modelingservices.testDescriptors.domain.Descriptors.DescriptorFactory.CheckAtomContainer;
import gov.epa.ccte.apps.modelingservices.testDescriptors.domain.Descriptors.DescriptorFactory.DescriptorData;
import gov.epa.ccte.apps.modelingservices.testDescriptors.domain.Descriptors.DescriptorFactory.DescriptorFactory;
import gov.epa.ccte.apps.modelingservices.testDescriptors.domain.Descriptors.DescriptorUtilities.HueckelAromaticityDetector;


import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.Vector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.ReaderFactory;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.springframework.stereotype.Service;

@Service
public class TestDescriptorService {

	private static final Logger logger = LogManager.getLogger(TestDescriptorService.class);

	public static String convertToTSV(Vector<String>vec) throws IOException {
		
		String result="";
		
		for (int i=0;i<vec.size();i++) {
			result+=vec.get(i);
			if(i<vec.size()-1) result+="\t";
		}
		return result;
	}


	static public IAtomContainer readOneMoleculeInAnyFormat(String anyFormatData) throws Exception {

		Reader r = null;
		ISimpleChemObjectReader reader = null;
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
//		try {
			r = new StringReader(anyFormatData);

			ReaderFactory readerFactory = new ReaderFactory();
			reader = readerFactory.createReader(r);

			IAtomContainer mol;

			if (reader == null){
				SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//				mol = smilesParser.parseSmiles(anyFormatData);
				
				try{
					mol = smilesParser.parseSmiles(anyFormatData);
				}catch (InvalidSmilesException ee){
//					System.out.println("here");
					throw ee;
				};
			}else
				mol = reader.read(builder.newAtomContainer());

			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
			try{
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
			}catch(Error e){
				System.out.println(e.getMessage());
			}
			StructureDiagramGenerator gen = new StructureDiagramGenerator();
			gen.generateCoordinates(mol);
			
			r.close();			
			if(reader!=null)reader.close();

			return mol;
			

//		}catch(Exception e){
//			throw e;
//		}
		
	}
	/**
	 * Generates 2d and 3d descriptors for molecule
	 * 
	 * @param sdf - molecule as MDL SDF with 3d coordinates included
	 * @return
	 * @throws IOException 
	 */
	private static Vector<String> goDescriptors2d_and_3dVector(String sdf, boolean also3D) throws Exception  {

		IAtomContainer ac=readOneMoleculeInAnyFormat(sdf); //I made local copy of this molecule reader 
		HueckelAromaticityDetector.debug=false;//suppress logging
		CheckAtomContainer.checkAtomContainer(ac, also3D);//need to check if can actually run descriptor calculations...

		DescriptorData dd=new DescriptorData(also3D);

		long t1 = System.currentTimeMillis();

		DescriptorFactory df = new DescriptorFactory(false);
		df.Calculate3DDescriptors=also3D;

		int descresult = df.CalculateDescriptors(ac, dd, true);

		logger.debug("Calculated descriptors in {}s", (System.currentTimeMillis() - t1) / 1000.);

		if (descresult == -1) 
			throw new Exception(df.errorMsg);

		return dd.getDescriptorValues();
	}

	public static Vector<String> goDescriptors2d_and_3dVector(String sdf) throws Exception {
		// TODO Auto-generated method stub
		return goDescriptors2d_and_3dVector(sdf,true);
	}

	public static Vector<String>  goDescriptors2dVector(String sdf) throws Exception {
		return goDescriptors2d_and_3dVector(sdf,false);
	}
	public static void main(String[] args) {
		try {
//			Vector<String>values=goDescriptors2dVector("CCCO");
//			Vector<String>values=goDescriptors2dVector("CCCOXX");
			Vector<String>values=goDescriptors2dVector("CB");
			System.out.println(values.get(0));
		} catch (Exception e) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
			System.out.println(e.getMessage().replace("\r", "").replace("\n", ""));
			
		}
	}

}
