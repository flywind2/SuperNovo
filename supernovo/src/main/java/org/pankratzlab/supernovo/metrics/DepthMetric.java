package org.pankratzlab.supernovo.metrics;

import java.util.function.Function;
import org.pankratzlab.supernovo.Pileup;
import org.pankratzlab.supernovo.Position;

public class DepthMetric extends SampleMetric {

  private enum DepthType {
    DEPTH("RawDepth", d -> d.getPileup(d.pos).get);

    private final String header;
    private final Function<DepthMetric, Number> depthCalc;

    private DepthType(String header, Function<DepthMetric, Number> depthCalc) {
      this.header = header;
      this.depthCalc = depthCalc;
    }
  }

  private final Position pos;
  private final byte a1;
  private final byte a2;

  /**
   * @param pileupFunc
   * @param pos
   * @param a1
   * @param a2
   */
  public DepthMetric(Function<Position, Pileup> pileupFunc, Position pos, byte a1, byte a2) {
    super(pileupFunc);
    this.pos = pos;
    this.a1 = a1;
    this.a2 = a2;
  }

  @Override
  public String getOutput() {
    // TODO Auto-generated method stub
    return null;
  }

  @Override
  public String getColumnHeader() {
    // TODO Auto-generated method stub
    return null;
  }
}
