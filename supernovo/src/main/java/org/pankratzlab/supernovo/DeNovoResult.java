package org.pankratzlab.supernovo;

import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import org.pankratzlab.supernovo.pileup.Depth;
import org.pankratzlab.supernovo.pileup.Pileup;
import com.google.common.collect.ImmutableList;

public class DeNovoResult {

  public enum Output {
    CHILD_ID("ID", r -> r.getChild().getId()),
    CHR("Chr", r -> r.getPos().getContig()),
    POS("Position", r -> r.getPos().getPosition()),
    A1("Allele_1", r -> (char) r.getA1()),
    A2("Allele_2", r -> (char) r.getA2()),
    BIALLELIC("Biallelic_Heterozygote", DeNovoResult::isBiallelic),
    DENOVO("DeNovo", DeNovoResult::isDeNovo),
    SUPERNOVO("SuperNovo", DeNovoResult::superNovo),
    HAP_CONCORDANCE("Haplotype_Concordance", DeNovoResult::meanConcordance),
    OVERLAP_DENOVOS("De_Novo_Variants_Overlapping_Reads", r -> r.getHapResults().getOtherDeNovos()),
    OVERLAP_TRIALLELICS(
        "Triallelic_Variants_Overlapping_Reads", r -> r.getHapResults().getOtherTriallelics()),
    CHILD_RAW_DEPTH("Raw_Depth", r -> r.getChild().getDepth().rawTotalDepth()),
    CHILD_A1_RAW_DEPTH(
        "Allele_1_Raw_Depth", r -> r.getChild().getDepth().allelicRawDepth(r.getA1())),
    CHILD_A2_RAW_DEPTH(
        "Allele_2_Raw_Depth", r -> r.getChild().getDepth().allelicRawDepth(r.getA2())),
    A1_CLIPPED_READS(
        "Allele_1_Clipped_Reads",
        r -> r.getChild().getPileup().getClippedReadCounts().count(r.getA1())),
    A2_CLIPPED_READS(
        "Allele_2_Clipped_Reads",
        r -> r.getChild().getPileup().getClippedReadCounts().count(r.getA2())),
    A1_UNMAPPED_MATE_READS(
        "Allele_1_Unmapped_Mate_Reads",
        r -> r.getChild().getPileup().getUnmappedMateCounts().count(r.getA1())),
    A2_UNMAPPED_MATE_READS(
        "Allele_2_Unmapped_Mate_Reads",
        r -> r.getChild().getPileup().getUnmappedMateCounts().count(r.getA2())),
    CHILD_WEIGHTED_DEPTH("Weighted_Depth", r -> r.getChild().getDepth().weightedTotalDepth()),
    CHILD_A1_WEIGHTED_DEPTH(
        "Allele_1_Weighted_Depth", r -> r.getChild().getDepth().allelicWeightedDepth(r.getA1())),
    CHILD_A2_WEIGHTED_DEPTH(
        "Allele_2_Weighted_Depth", r -> r.getChild().getDepth().allelicWeightedDepth(r.getA2())),
    P1_ID("Parent1_ID", r -> r.getParent1().getId()),
    P1_RAW_DEPTH("Parent1_Raw_Depth", r -> r.getParent1().getDepth().rawTotalDepth()),
    P1_A1_RAW_DEPTH(
        "Parent1_Allele_1_Raw_Depth", r -> r.getParent1().getDepth().allelicRawDepth(r.getA1())),
    P1_A2_RAW_DEPTH(
        "Parent1_Allele_2_Raw_Depth", r -> r.getParent1().getDepth().allelicRawDepth(r.getA2())),
    P1_WEIGHTED_DEPTH(
        "Parent1_Weighted_Depth", r -> r.getParent1().getDepth().weightedTotalDepth()),
    P1_A1_WEIGHTED_DEPTH(
        "Parent1_Allele_1_Weighted_Depth",
        r -> r.getParent1().getDepth().allelicWeightedDepth(r.getA1())),
    P1_A2_WEIGHTED_DEPTH(
        "Parent1_Allele_2_Weighted_Depth",
        r -> r.getParent1().getDepth().allelicWeightedDepth(r.getA2())),
    P2_ID("Parent2_ID", r -> r.getParent2().getId()),
    P2_RAW_DEPTH("Parent2_Raw_Depth", r -> r.getParent2().getDepth().rawTotalDepth()),
    P2_A1_RAW_DEPTH(
        "Parent2_Allele_1_Raw_Depth", r -> r.getParent2().getDepth().allelicRawDepth(r.getA1())),
    P2_A2_RAW_DEPTH(
        "Parent2_Allele_2_Raw_Depth", r -> r.getParent2().getDepth().allelicRawDepth(r.getA2())),
    P2_WEIGHTED_DEPTH(
        "Parent2_Weighted_Depth", r -> r.getParent2().getDepth().weightedTotalDepth()),
    P2_A1_WEIGHTED_DEPTH(
        "Parent2_Allele_1_Weighted_Depth",
        r -> r.getParent2().getDepth().allelicWeightedDepth(r.getA1())),
    P2_A2_WEIGHTED_DEPTH(
        "Parent2_Allele_2_Weighted_Depth",
        r -> r.getParent2().getDepth().allelicWeightedDepth(r.getA2()));

    private static final String DELIM = "\t";

    private final String columnHeader;
    private final Function<DeNovoResult, String> dataGenerator;

    Output(String columnHeader, Function<DeNovoResult, Object> getter) {
      this.columnHeader = columnHeader;
      dataGenerator = getter.andThen(Object::toString);
    }

    /** @return the columnHeader */
    public String getColumnHeader() {
      return columnHeader;
    }

    public String getData(DeNovoResult result) {
      return dataGenerator.apply(result);
    }

    public static String generateHeaderLine() {
      return Arrays.stream(values())
          .map(Output::getColumnHeader)
          .collect(Collectors.joining(DELIM));
    }

    public static String generateOutputLine(DeNovoResult result) {
      return Arrays.stream(values()).map(o -> o.getData(result)).collect(Collectors.joining(DELIM));
    }
  }

  public static class Sample {

    private final String id;
    private final Pileup pileup;

    /**
     * @param id
     * @param pileup
     */
    public Sample(String id, Pileup pileup) {
      super();
      this.id = id;
      this.pileup = pileup;
    }

    /** @return the id */
    public String getId() {
      return id;
    }

    /** @return the depth */
    public Depth getDepth() {
      return pileup.getDepth();
    }

    /** @return the pileup */
    public Pileup getPileup() {
      return pileup;
    }
  }

  private final Position pos;
  private final HaplotypeEvaluator.Result hapResults;
  private final boolean biallelic;
  private final boolean deNovo;
  private final Sample child;
  private final List<Sample> parents;

  public DeNovoResult(
      Position pos, HaplotypeEvaluator.Result hapResults, Sample child, Sample p1, Sample p2) {
    this.pos = pos;
    this.hapResults = hapResults;
    this.child = child;
    this.parents = ImmutableList.of(p1, p2);
    this.biallelic = TrioEvaluator.looksBiallelic(child.getPileup());
    this.deNovo = TrioEvaluator.looksDenovo(child.getPileup(), p1.getPileup(), p2.getPileup());
  }

  /** @return the pos */
  public Position getPos() {
    return pos;
  }

  /** @return the a1 */
  public byte getA1() {
    return child.getDepth().getA1().orElseThrow(IllegalStateException::new);
  }

  /** @return the a2 */
  public byte getA2() {
    return child.getDepth().getA2().orElseThrow(IllegalStateException::new);
  }

  /** @return the child */
  public Sample getChild() {
    return child;
  }

  /** @return the parents */
  public List<Sample> getParents() {
    return parents;
  }

  /** @return the parent1 */
  public Sample getParent1() {
    return parents.get(0);
  }

  /** @return the parent2 */
  public Sample getParent2() {
    return parents.get(1);
  }

  /** @return the hapResults */
  public HaplotypeEvaluator.Result getHapResults() {
    return hapResults;
  }

  public double meanConcordance() {
    if (getHapResults().getConcordances().isEmpty()) return 1.0;
    return getHapResults()
        .getConcordances()
        .stream()
        .mapToDouble(Double::valueOf)
        .summaryStatistics()
        .getAverage();
  }

  /** @return the biallelic */
  public boolean isBiallelic() {
    return biallelic;
  }

  /** @return the deNovo */
  public boolean isDeNovo() {
    return deNovo;
  }

  public boolean superNovo() {
    return biallelic
        && deNovo
        && hapResults.getOtherDeNovos() == 0
        && meanConcordance() >= 0.95
        && hapResults.getOtherTriallelics() == 0;
  }
}
