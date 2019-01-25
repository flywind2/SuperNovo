package org.pankratzlab.supernovo.output;

import java.util.List;
import java.util.Optional;
import org.pankratzlab.supernovo.HaplotypeEvaluator;
import org.pankratzlab.supernovo.Position;
import org.pankratzlab.supernovo.TrioEvaluator;
import org.pankratzlab.supernovo.HaplotypeEvaluator.Result;
import org.pankratzlab.supernovo.pileup.Depth;
import org.pankratzlab.supernovo.pileup.Pileup;
import com.google.common.collect.ImmutableList;

public class DeNovoResult implements OutputFields {

  public static class Sample implements OutputFields {

    public final String id;
    public final int rawDepth;
    public final int a1RawDepth;
    public final int a2RawDepth;
    public final int a1ClippedReads;
    public final int a2ClippedReads;
    public final int a1UnmappedMateReads;
    public final int a2UnmappedMateReads;
    public final double weightedDepth;
    public final double a1WeightedDepth;
    public final double a2WeightedDepth;

    private final Pileup pileup;

    /**
     * @param id
     * @param pileup
     */
    public Sample(String id, Pileup pileup, Optional<Byte> a1, Optional<Byte> a2) {
      super();
      this.pileup = pileup;
      Depth depth = pileup.getDepth();

      this.id = id;
      rawDepth = depth.rawTotalDepth();
      a1RawDepth = a1.map(depth::allelicRawDepth).orElse(0);
      a2RawDepth = a2.map(depth::allelicRawDepth).orElse(0);
      a1ClippedReads = a1.map(pileup.getClippedReadCounts()::count).orElse(0);
      a2ClippedReads = a2.map(pileup.getClippedReadCounts()::count).orElse(0);
      a1UnmappedMateReads = a1.map(pileup.getUnmappedMateCounts()::count).orElse(0);
      a2UnmappedMateReads = a2.map(pileup.getUnmappedMateCounts()::count).orElse(0);
      weightedDepth = depth.weightedTotalDepth();
      a1WeightedDepth = a1.map(depth::allelicWeightedDepth).orElse(0.0);
      a2WeightedDepth = a2.map(depth::allelicWeightedDepth).orElse(0.0);
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

  public final String chr;
  public final int position;
  public final char allele1;
  public final char allele2;
  public final boolean biallelicHeterozygote;
  public final boolean deNovo;
  public final boolean superNovo;
  public final double meanHaplotypeConcordance;
  public final int overlappingReadsDeNovoCount;
  public final int overlapingReadsThirdAlleleCount;
  public final Sample child;
  public final Sample p1;
  public final Sample p2;

  private final Position pos;
  private final HaplotypeEvaluator.Result hapResults;
  private final List<Sample> parents;

  public DeNovoResult(
      Position pos, HaplotypeEvaluator.Result hapResults, Sample child, Sample p1, Sample p2) {
    this.pos = pos;
    this.hapResults = hapResults;
    this.child = child;
    this.p1 = p1;
    this.p2 = p2;
    this.parents = ImmutableList.of(p1, p2);

    position = pos.getPosition();
    chr = pos.getContig();
    allele1 = child.getDepth().getA1().map(b -> (char) b.byteValue()).orElse('.');
    allele2 = child.getDepth().getA2().map(b -> (char) b.byteValue()).orElse('.');
    biallelicHeterozygote = TrioEvaluator.looksBiallelic(child.getPileup());
    deNovo = TrioEvaluator.looksDenovo(child.getPileup(), p1.getPileup(), p2.getPileup());
    if (hapResults.getConcordances().isEmpty()) meanHaplotypeConcordance = 1.0;
    else
      meanHaplotypeConcordance =
          hapResults
              .getConcordances()
              .stream()
              .mapToDouble(Double::valueOf)
              .summaryStatistics()
              .getAverage();
    superNovo =
        biallelicHeterozygote
            && deNovo
            && hapResults.getOtherDeNovos() == 0
            && meanHaplotypeConcordance >= 0.95
            && hapResults.getOtherTriallelics() == 0;
    overlappingReadsDeNovoCount = hapResults.getOtherDeNovos();
    overlapingReadsThirdAlleleCount = hapResults.getOtherTriallelics();
  }
}
