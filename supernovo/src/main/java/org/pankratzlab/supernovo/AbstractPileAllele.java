package org.pankratzlab.supernovo;

public abstract class AbstractPileAllele implements PileAllele {

  private final String alleleString;

  /** @param alleleString */
  public AbstractPileAllele(String alleleString) {
    super();
    this.alleleString = alleleString;
  }

  @Override
  public String toString() {
    return alleleString;
  }
}
