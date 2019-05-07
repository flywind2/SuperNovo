package org.pankratzlab.supernovo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import com.google.common.collect.ImmutableRangeSet;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.vcf.VCFFileReader;
import picocli.CommandLine;
import picocli.CommandLine.ArgGroup;
import picocli.CommandLine.Option;
import picocli.CommandLine.RunAll;

public class App implements Runnable {

  private static class Query {
    @Option(
      names = {"--vcf", "-v"},
      paramLabel = "VCF",
      description = "VCF with variants to query for de novo mutations",
      required = true
    )
    private File vcf;

    @Option(
      names = {"--fasta", "-f"},
      paramLabel = "FASTA",
      description = "Reference genome FASTA to query for de novo mutations",
      required = true
    )
    private File referenceFasta;
  }

  @ArgGroup(exclusive = true, multiplicity = "1")
  private Query query;

  @Option(
    names = {"--childBam", "--bam"},
    paramLabel = "BAM",
    description = "BAM of child",
    required = true
  )
  private File childBam;

  @Option(
    names = {"--childID", "--cID"},
    paramLabel = "ID",
    description = "Sample ID of child",
    required = true
  )
  private String childID;

  @Option(
    names = {"--parent1Bam", "--p1Bam"},
    paramLabel = "BAM",
    description = "BAM of parent 1",
    required = true
  )
  private File p1Bam;

  @Option(
    names = {"--parent1ID", "--p1ID"},
    paramLabel = "ID",
    description = "Sample ID of parent 1",
    required = true
  )
  private String p1ID;

  @Option(
    names = {"--parent2Bam", "--p2Bam"},
    paramLabel = "BAM",
    description = "BAM of parent 2",
    required = true
  )
  private File p2Bam;

  @Option(
    names = {"--parent2ID", "--p2ID"},
    paramLabel = "ID",
    description = "Sample ID of parent 2",
    required = true
  )
  private String p2ID;

  @Option(
    names = {"--white-list", "-i"},
    paramLabel = "FILE",
    description = "White list region file"
  )
  private RangeSet<GenomePosition> whiteList = ImmutableRangeSet.of(Range.all());

  @Option(
    names = {"--black-list", "-x"},
    paramLabel = "FILE",
    description = "Black list region file"
  )
  private RangeSet<GenomePosition> blackList = ImmutableRangeSet.of();

  @Option(
    names = {"--output", "-o"},
    paramLabel = "FILE",
    description = "Output file for parsed de novo variants",
    required = true
  )
  private File output;

  public static void main(String[] args) {
    CommandLine cl = new CommandLine(new App());
    cl.registerConverter(RangeSet.class, s -> bedToRangeSet(new File(s)));
    cl.parseWithHandler(new RunAll(), args);
  }

  private static RangeSet<GenomePosition> bedToRangeSet(File bedFile)
      throws FileNotFoundException, IOException {
    BEDCodec bedCodec = new BEDCodec();
    try (BufferedReader bedReader = new BufferedReader(new FileReader(bedFile))) {
      return bedReader
          .lines()
          .map(bedCodec::decodeLoc)
          .map(App::locToRange)
          .collect(TreeRangeSet::create, TreeRangeSet::add, TreeRangeSet::addAll);
    }
  }

  private static Range<GenomePosition> locToRange(BEDFeature loc) {
    return Range.closed(
        new GenomePosition(loc.getContig(), loc.getStart()),
        new GenomePosition(loc.getContig(), loc.getEnd()));
  }

  private RangeSet<GenomePosition> generateTargetIntervals() {
    RangeSet<GenomePosition> targetIntervals = TreeRangeSet.create(whiteList);
    targetIntervals.removeAll(blackList);
    return ImmutableRangeSet.copyOf(targetIntervals);
  }

  @Override
  public void run() {
    SamReaderFactory srFactory = SamReaderFactory.make();
    try (SamReader child = srFactory.open(childBam);
        SamReader childReader2 = srFactory.open(childBam);
        SamReader p1 = srFactory.open(p1Bam);
        SamReader p2 = srFactory.open(p2Bam);
        TrioEvaluator trioEvaluator =
            new TrioEvaluator(
                child, childReader2, childID, p1, p1ID, p2, p2ID, generateTargetIntervals())) {
      if (query.vcf != null) {
        try (VCFFileReader vcfReader = new VCFFileReader(query.vcf)) {
          trioEvaluator.reportDeNovos(vcfReader, output);
        }
      } else {
        try (IndexedFastaSequenceFile fastaSeq =
            new IndexedFastaSequenceFile(query.referenceFasta)) {
          trioEvaluator.reportAllDeNovos(fastaSeq, output);
        }
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
