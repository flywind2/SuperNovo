package org.pankratzlab.supernovo;

import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import com.google.common.collect.ImmutableList;

public class DeNovoResult {

  public enum Output {
    A1("Allele_1", DeNovoResult::getA1),
    A2("Allele_2", DeNovoResult::getA2),
    CHILD_DEPTH("Depth", r -> r.getChild().getDepth()),
    CHILD_A1_DEPTH("Allele_1_Depth", r -> r.getChild().getA1Depth()),
    CHILD_A2_DEPTH("Allele_2_Depth", r -> r.getChild().getA2Depth()),
    P1_DEPTH("Parent1_Depth", r -> r.getParent1().getDepth()),
    P1_A1_DEPTH("Parent1_Allele_1_Depth", r -> r.getParent1().getA1Depth()),
    P1_A2_DEPTH("Parent1_Allele_2_Depth", r -> r.getParent1().getA2Depth()),
    P2_DEPTH("Parent2_Depth", r -> r.getParent2().getDepth()),
    P2_A1_DEPTH("Parent2_Allele_1_Depth", r -> r.getParent2().getA1Depth()),
    P2_A2_DEPTH("Parent2_Allele_2_Depth", r -> r.getParent2().getA2Depth());

    private static final String DELIM = "\t";

    private final String columnHeader;
    private final Function<DeNovoResult, String> dataGenerator;

    Output(String columnHeader, Function<DeNovoResult, Object> getter) {
      this.columnHeader = columnHeader;
      dataGenerator = getter.andThen(Object::toString);
    }

    /**
     * @return the columnHeader
     */
    public String getColumnHeader() {
      return columnHeader;
    }

    public String getData(DeNovoResult result) {
      return dataGenerator.apply(result);
    }

    public static String generateHeaderLine() {
      return Arrays.stream(values()).map(Output::getColumnHeader)
                   .collect(Collectors.joining(DELIM));
    }

    public static String generateOutputLine(DeNovoResult result) {
      return Arrays.stream(values()).map(o -> o.getData(result)).collect(Collectors.joining(DELIM));
    }

  }

  public static class Sample {

    private final String id;
    private final int depth;
    private final int a1Depth;
    private final int a2Depth;

    /**
     * @param id TODO
     * @param depth
     * @param a1Depth
     * @param a2Depth
     */
    public Sample(String id, int depth, int a1Depth, int a2Depth) {
      super();
      this.id = id;
      this.depth = depth;
      this.a1Depth = a1Depth;
      this.a2Depth = a2Depth;
    }

    /**
     * @return the depth
     */
    public int getDepth() {
      return depth;
    }

    /**
     * @return the a1Depth
     */
    public int getA1Depth() {
      return a1Depth;
    }

    /**
     * @return the a2Depth
     */
    public int getA2Depth() {
      return a2Depth;
    }

  }

  private final byte a1;
  private final byte a2;
  private final Sample child;
  private final List<Sample> parents;

  public DeNovoResult(byte a1, byte a2, Sample child, Sample p1, Sample p2) {
    this.a1 = a1;
    this.a2 = a2;
    this.child = child;
    this.parents = ImmutableList.of(p1, p2);
  }

  /**
   * @return the a1
   */
  public byte getA1() {
    return a1;
  }

  /**
   * @return the a2
   */
  public byte getA2() {
    return a2;
  }

  /**
   * @return the child
   */
  public Sample getChild() {
    return child;
  }

  /**
   * @return the parents
   */
  public List<Sample> getParents() {
    return parents;
  }

  /**
   * @return the parent1
   */
  public Sample getParent1() {
    return parents.get(0);
  }

  /**
   * @return the parent2
   */
  public Sample getParent2() {
    return parents.get(1);
  }

}
