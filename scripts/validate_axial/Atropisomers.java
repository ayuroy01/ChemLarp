import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.stereo.*;
import org.openscience.cdk.CDKConstants;
import java.io.*;
public class Atropisomers { public static void main(String[] a) throws Exception {
  IteratingSDFReader r = new IteratingSDFReader(new FileReader(a[0]), SilentChemObjectBuilder.getInstance());
  int n=0, ok=0, noaxis=0;
  while (r.hasNext()) { IAtomContainer m = r.next(); n++;
    String ours = m.getTitle().split(" ")[0];
    m.setStereoElements(StereoElementFactory.using3DCoordinates(m).createAll());
    com.simolecule.centres.CdkLabeller.label(m);
    String lab = null;
    for (IStereoElement<?,?> se : m.stereoElements()) if (se instanceof Atropisomeric) { IBond b = (IBond) se.getFocus(); lab = b.getProperty(CDKConstants.CIP_DESCRIPTOR) == null ? "none" : b.getProperty(CDKConstants.CIP_DESCRIPTOR).toString(); }
    if (lab == null) lab = "none"; if (ours.equals("none") && lab.equals("none")) { noaxis++; continue; }
    String mapped = lab.equals("M") ? "Ra" : lab.equals("P") ? "Sa" : lab;
    if (mapped.equals(ours)) ok++; else System.out.println("DIFF " + m.getTitle() + " cdk " + lab);
  }
  System.out.printf("%d biaryls: %d with the same Ra/Sa, %d not stereogenic in both, %d different%n", n, ok, noaxis, n - ok - noaxis); } }
