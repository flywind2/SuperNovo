package org.pankratzlab.supernovo.utilities;

public final class Phred {

  private Phred() {}

  public static double getErrorProbability(int phred) {
    return Math.pow(10, Math.negateExact(phred) / 10.0);
  }

  public static double getAccuracy(int phred) {
    return 1.0 - getErrorProbability(phred);
  }
}
