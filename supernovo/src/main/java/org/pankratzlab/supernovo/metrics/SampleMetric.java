package org.pankratzlab.supernovo.metrics;

public abstract class SampleMetric implements Metric {
  public enum Type {
    CHILD,
    P1,
    P2;
  }

  private final Type type;
  private final String baseHeader;

  /**
   * @param type
   * @param baseHeader
   */
  public SampleMetric(Type type, String baseHeader) {
    super();
    this.type = type;
    this.baseHeader = baseHeader;
  }

  @Override
  public String getColumnHeader() {
    return type.toString() + baseHeader;
  }
}
