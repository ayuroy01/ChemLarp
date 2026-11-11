// PDB fixtures (public domain, RCSB): 1CRN residues 1-5 (crambin, heavy atoms),
// 1L2Y model 1 residues 1-3 (Trp-cage NMR, with hydrogens) plus an empty second model.
const char* const kCrambin5 = R"PDB(HEADER    PLANT SEED PROTEIN                      30-APR-81   1CRN
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N  
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C  
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C  
ATOM      4  O   THR A   1      15.268  13.825   5.594  1.00  9.85           O  
ATOM      5  CB  THR A   1      18.170  12.703   5.337  1.00 13.02           C  
ATOM      6  OG1 THR A   1      19.334  12.829   4.463  1.00 15.06           O  
ATOM      7  CG2 THR A   1      18.150  11.546   6.304  1.00 14.23           C  
ATOM      8  N   THR A   2      15.115  11.555   5.265  1.00  7.81           N  
ATOM      9  CA  THR A   2      13.856  11.469   6.066  1.00  8.31           C  
ATOM     10  C   THR A   2      14.164  10.785   7.379  1.00  5.80           C  
ATOM     11  O   THR A   2      14.993   9.862   7.443  1.00  6.94           O  
ATOM     12  CB  THR A   2      12.732  10.711   5.261  1.00 10.32           C  
ATOM     13  OG1 THR A   2      13.308   9.439   4.926  1.00 12.81           O  
ATOM     14  CG2 THR A   2      12.484  11.442   3.895  1.00 11.90           C  
ATOM     15  N   CYS A   3      13.488  11.241   8.417  1.00  5.24           N  
ATOM     16  CA  CYS A   3      13.660  10.707   9.787  1.00  5.39           C  
ATOM     17  C   CYS A   3      12.269  10.431  10.323  1.00  4.45           C  
ATOM     18  O   CYS A   3      11.393  11.308  10.185  1.00  6.54           O  
ATOM     19  CB  CYS A   3      14.368  11.748  10.691  1.00  5.99           C  
ATOM     20  SG  CYS A   3      15.885  12.426  10.016  1.00  7.01           S  
ATOM     21  N   CYS A   4      12.019   9.272  10.928  1.00  3.90           N  
ATOM     22  CA  CYS A   4      10.646   8.991  11.408  1.00  4.24           C  
ATOM     23  C   CYS A   4      10.654   8.793  12.919  1.00  3.72           C  
ATOM     24  O   CYS A   4      11.659   8.296  13.491  1.00  5.30           O  
ATOM     25  CB  CYS A   4      10.057   7.752  10.682  1.00  4.41           C  
ATOM     26  SG  CYS A   4       9.837   8.018   8.904  1.00  4.72           S  
ATOM     27  N   PRO A   5       9.561   9.108  13.563  1.00  3.96           N  
ATOM     28  CA  PRO A   5       9.448   9.034  15.012  1.00  4.25           C  
ATOM     29  C   PRO A   5       9.288   7.670  15.606  1.00  4.96           C  
ATOM     30  O   PRO A   5       9.490   7.519  16.819  1.00  7.44           O  
ATOM     31  CB  PRO A   5       8.230   9.957  15.345  1.00  5.11           C  
ATOM     32  CG  PRO A   5       7.338   9.786  14.114  1.00  5.24           C  
ATOM     33  CD  PRO A   5       8.366   9.804  12.958  1.00  5.20           C  
END
)PDB";
const char* const kTrpCage3 = R"PDB(MODEL        1
ATOM      1  N   ASN A   1      -8.901   4.127  -0.555  1.00  0.00           N  
ATOM      2  CA  ASN A   1      -8.608   3.135  -1.618  1.00  0.00           C  
ATOM      3  C   ASN A   1      -7.117   2.964  -1.897  1.00  0.00           C  
ATOM      4  O   ASN A   1      -6.634   1.849  -1.758  1.00  0.00           O  
ATOM      5  CB  ASN A   1      -9.437   3.396  -2.889  1.00  0.00           C  
ATOM      6  CG  ASN A   1     -10.915   3.130  -2.611  1.00  0.00           C  
ATOM      7  OD1 ASN A   1     -11.269   2.700  -1.524  1.00  0.00           O  
ATOM      8  ND2 ASN A   1     -11.806   3.406  -3.543  1.00  0.00           N  
ATOM      9  H1  ASN A   1      -8.330   3.957   0.261  1.00  0.00           H  
ATOM     10  H2  ASN A   1      -8.740   5.068  -0.889  1.00  0.00           H  
ATOM     11  H3  ASN A   1      -9.877   4.041  -0.293  1.00  0.00           H  
ATOM     12  HA  ASN A   1      -8.930   2.162  -1.239  1.00  0.00           H  
ATOM     13  HB2 ASN A   1      -9.310   4.417  -3.193  1.00  0.00           H  
ATOM     14  HB3 ASN A   1      -9.108   2.719  -3.679  1.00  0.00           H  
ATOM     15 HD21 ASN A   1     -11.572   3.791  -4.444  1.00  0.00           H  
ATOM     16 HD22 ASN A   1     -12.757   3.183  -3.294  1.00  0.00           H  
ATOM     17  N   LEU A   2      -6.379   4.031  -2.228  1.00  0.00           N  
ATOM     18  CA  LEU A   2      -4.923   4.002  -2.452  1.00  0.00           C  
ATOM     19  C   LEU A   2      -4.136   3.187  -1.404  1.00  0.00           C  
ATOM     20  O   LEU A   2      -3.391   2.274  -1.760  1.00  0.00           O  
ATOM     21  CB  LEU A   2      -4.411   5.450  -2.619  1.00  0.00           C  
ATOM     22  CG  LEU A   2      -4.795   6.450  -1.495  1.00  0.00           C  
ATOM     23  CD1 LEU A   2      -3.612   6.803  -0.599  1.00  0.00           C  
ATOM     24  CD2 LEU A   2      -5.351   7.748  -2.084  1.00  0.00           C  
ATOM     25  H   LEU A   2      -6.821   4.923  -2.394  1.00  0.00           H  
ATOM     26  HA  LEU A   2      -4.750   3.494  -3.403  1.00  0.00           H  
ATOM     27  HB2 LEU A   2      -3.340   5.414  -2.672  1.00  0.00           H  
ATOM     28  HB3 LEU A   2      -4.813   5.817  -3.564  1.00  0.00           H  
ATOM     29  HG  LEU A   2      -5.568   6.022  -0.858  1.00  0.00           H  
ATOM     30 HD11 LEU A   2      -3.207   5.905  -0.146  1.00  0.00           H  
ATOM     31 HD12 LEU A   2      -2.841   7.304  -1.183  1.00  0.00           H  
ATOM     32 HD13 LEU A   2      -3.929   7.477   0.197  1.00  0.00           H  
ATOM     33 HD21 LEU A   2      -4.607   8.209  -2.736  1.00  0.00           H  
ATOM     34 HD22 LEU A   2      -6.255   7.544  -2.657  1.00  0.00           H  
ATOM     35 HD23 LEU A   2      -5.592   8.445  -1.281  1.00  0.00           H  
ATOM     36  N   TYR A   3      -4.354   3.455  -0.111  1.00  0.00           N  
ATOM     37  CA  TYR A   3      -3.690   2.738   0.981  1.00  0.00           C  
ATOM     38  C   TYR A   3      -4.102   1.256   1.074  1.00  0.00           C  
ATOM     39  O   TYR A   3      -3.291   0.409   1.442  1.00  0.00           O  
ATOM     40  CB  TYR A   3      -3.964   3.472   2.302  1.00  0.00           C  
ATOM     41  CG  TYR A   3      -2.824   3.339   3.290  1.00  0.00           C  
ATOM     42  CD1 TYR A   3      -2.746   2.217   4.138  1.00  0.00           C  
ATOM     43  CD2 TYR A   3      -1.820   4.326   3.332  1.00  0.00           C  
ATOM     44  CE1 TYR A   3      -1.657   2.076   5.018  1.00  0.00           C  
ATOM     45  CE2 TYR A   3      -0.725   4.185   4.205  1.00  0.00           C  
ATOM     46  CZ  TYR A   3      -0.639   3.053   5.043  1.00  0.00           C  
ATOM     47  OH  TYR A   3       0.433   2.881   5.861  1.00  0.00           O  
ATOM     48  H   TYR A   3      -4.934   4.245   0.120  1.00  0.00           H  
ATOM     49  HA  TYR A   3      -2.615   2.768   0.796  1.00  0.00           H  
ATOM     50  HB2 TYR A   3      -4.117   4.513   2.091  1.00  0.00           H  
ATOM     51  HB3 TYR A   3      -4.886   3.096   2.750  1.00  0.00           H  
ATOM     52  HD1 TYR A   3      -3.513   1.456   4.101  1.00  0.00           H  
ATOM     53  HD2 TYR A   3      -1.877   5.200   2.695  1.00  0.00           H  
ATOM     54  HE1 TYR A   3      -1.576   1.221   5.669  1.00  0.00           H  
ATOM     55  HE2 TYR A   3       0.033   4.952   4.233  1.00  0.00           H  
ATOM     56  HH  TYR A   3       1.187   3.395   5.567  1.00  0.00           H  
ENDMDL
MODEL        2
ENDMDL
END
)PDB";
