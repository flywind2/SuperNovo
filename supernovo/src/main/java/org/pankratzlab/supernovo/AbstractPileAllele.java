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

  /* (non-Javadoc)
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((alleleString == null) ? 0 : alleleString.hashCode());
    return result;
  }

  /* (non-Javadoc)
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj instanceof AbstractPileAllele)) return false;
    AbstractPileAllele other = (AbstractPileAllele) obj;
    if (alleleString == null) {
      if (other.alleleString != null) return false;
    } else if (!alleleString.equals(other.alleleString)) return false;
    return true;
  }
}
