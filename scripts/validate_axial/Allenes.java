import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.*;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.stereo.*;
import java.nio.file.*;
import java.util.*;

public class Allenes {
  static SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
  static SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Absolute);
  static String label(IAtomContainer m) {
    com.simolecule.centres.CdkLabeller.label(m);
    StringBuilder b = new StringBuilder();
    for (IStereoElement<?,?> se : m.stereoElements()) {
      if (se instanceof ExtendedTetrahedral || se instanceof Atropisomeric) {
        Object f = se.getFocus();
        Object d = (f instanceof IAtom) ? ((IAtom) f).getProperty(CDKConstants.CIP_DESCRIPTOR) : ((IBond) f).getProperty(CDKConstants.CIP_DESCRIPTOR);
        b.append(d);
      }
    }
    return b.toString();
  }
  public static void main(String[] a) throws Exception {
    int n = 0, labOk = 0, canOk = 0, rndOk = 0, rnd = 0;
    Map<String,String> map = new HashMap<>();
    for (String line : Files.readAllLines(Path.of(a[0]))) {
      String[] t = line.split("\t");
      IAtomContainer m = sp.parseSmiles(t[0]);
      String ref = sg.create(m);
      String cl = label(m);
      String ours = t[1];
      n++;
      String mapped = cl.equals("R") ? "Ra" : cl.equals("S") ? "Sa" : cl.equals("M") ? "Ra" : cl.equals("P") ? "Sa" : cl;
      map.put(cl, ours);
      if (mapped.equals(ours)) labOk++; else System.out.println("LABEL " + t[0] + " ours " + ours + " cdk " + cl);
      if (sg.create(sp.parseSmiles(t[2])).equals(ref)) canOk++; else System.out.println("CANON " + t[0] + " -> " + t[2]);
      for (int k = 3; k < t.length; k++) { rnd++; if (sg.create(sp.parseSmiles(t[k])).equals(ref)) rndOk++; else System.out.println("RANDOM " + t[0] + " -> " + t[k]); }
    }
    System.out.println("cdk label symbols seen: " + map.keySet());
    System.out.printf("%d allenes: label %d, canonical SMILES same stereo per CDK %d, random spellings %d/%d%n", n, labOk, canOk, rndOk, rnd);
  }
}
