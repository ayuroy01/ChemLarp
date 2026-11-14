import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.*;
import java.nio.file.*;
public class WriteCdk { public static void main(String[] a) throws Exception {
  SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
  SmilesGenerator abs = new SmilesGenerator(SmiFlavor.Absolute), stereo = new SmilesGenerator(SmiFlavor.Isomeric);
  for (String line : Files.readAllLines(Path.of(a[0]))) { String[] t = line.split("\t"); IAtomContainer m = sp.parseSmiles(t[0]);
    System.out.println(t[0] + "\t" + t[1] + "\t" + abs.create(m) + "\t" + stereo.create(m)); } } }
